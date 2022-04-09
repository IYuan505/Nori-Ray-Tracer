#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/vector.h>

NORI_NAMESPACE_BEGIN

class PathMisIntegrator : public Integrator {
public:
    PathMisIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;           // Intersection information
        Intersection itsBRDF;       // Intersection information for BRDF MIS
        Ray3f iterRay(ray);         // Ray used for path tracing
        Color3f throughput = {1};   // Path throughput
        Color3f fr;                 // Reflection or refraction factor (albedo)
        Color3f L = {0};            // Final rendered light 
        Color3f light;              // MIS that come from direct emitter sampling
        Color3f brdf;               // MIS that come from BRDF  
        Color3f emitterLight = {0}; // Light of emitter
        Vector3f pToLight;          // Vector from point to light       
        std::vector<Mesh *> meshes; // All meshes of the scene
        EmitterQueryRecord eQ;      // The emitter query record
        float eta = 1.0f;           // Cumulative eta
        float prob;                 // Probability to continue
        float lightProb;            // Probability of sampling on light
        float lightPortion;         // Light portion from MIS
        float brdfProb;             // Probability of sampling the BRDF direction
        float brdfPortion;          // BRDF portion from MIS
        bool foundIntersection;     // Whether there is an intersection point found
        bool specularBounce = false;// Whether last bounce is a specular bounce

        /* Start the path tracing */
        for (uint32_t bounces = 0; ; ++bounces) {
            /* If there is no intersection found, break from the loop */
            foundIntersection = scene->rayIntersect(iterRay, its);
            if (!foundIntersection) break;

            /* Add the emitter light at path vertex */
            if (its.mesh->isEmitter()) {
                /* Avoid double counting but also consider MIS
                   bQ comes from the previous iteration */
                lightProb = 1 / its.mesh->getArea() * its.t * its.t / -iterRay.d.dot(its.shFrame.n);
                brdfPortion = (bounces == 0 || specularBounce) ? 1 : brdfProb / (lightProb + brdfProb);
                eQ = EmitterQueryRecord(-iterRay.d, its.shFrame.n);
                L += throughput * its.mesh->getEmitter()->eval(eQ) * brdfPortion;
            }
            
            /* Multiple importance sampling (MIS):
               Sampling direct emitters, rectified by d^2/theta */
            light = scene->uniformlySampleLight(sampler, &its, &eQ, iterRay);
            if (light.getLuminance() > 0) {
                /* Only consider visible light, cos > 0
                   Convert PDF in area to PDF in solid angles */
                lightProb = eQ.pdf * eQ.d * eQ.d / eQ.wo.dot(eQ.n);
                BSDFQueryRecord bQLight(its.shFrame.toLocal(-iterRay.d).normalized(), 
                    its.shFrame.toLocal(-eQ.wo).normalized(), ESolidAngle);
                brdfProb = its.mesh->getBSDF()->pdf(bQLight);
                lightPortion = lightProb / (lightProb + brdfProb);
                L += throughput * light * lightPortion;
            }

            /* Account for indirect light, we sample a new direction on this surface */
            BSDFQueryRecord bQ(its.shFrame.toLocal(-iterRay.d).normalized());
            fr = its.mesh->getBSDF()->sample(bQ, sampler->next2D());
            brdfProb = its.mesh->getBSDF()->pdf(bQ);
            if (fr.getLuminance() == 0.0f) break;

            /* Update the iteration ray given the sampling */
            Ray3f temp = Ray3f(its.p, its.shFrame.toWorld(bQ.wo));
            memcpy(&iterRay, &temp, sizeof(Ray3f));

            /* Compute the probability to continue */
            eta *= bQ.eta;
            prob = std::min(0.99f, throughput.maxCoeff() * eta * eta);

            /* Only start doing Russian Roulette after at least three bounces */
            prob = bounces < 4 ? 1 : prob;
            if (sampler->next1D() < prob) {
                throughput *= fr / prob;
                specularBounce = !its.mesh->getBSDF()->isDiffuse();
            } else {
                break;
            }
        }
        return L;
    }

    std::string toString() const {
        return "PathMisIntegrator[]";
    }

};

NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");
NORI_NAMESPACE_END