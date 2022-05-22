#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/vector.h>

NORI_NAMESPACE_BEGIN

class PathMatsIntegrator : public Integrator {
public:
    PathMatsIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;           // Intersection information
        Ray3f iterRay(ray);         // Ray used for path tracing
        Color3f throughput = {1};   // Path throughput
        Color3f fr;                 // Reflection or refraction factor (albedo)
        Color3f L = {0};            // Final rendered light 
        const BSDF *bsdf;           // BSDF of the current mesh
        EmitterQueryRecord eQ;      // The emitter query record
        float eta = 1.0f;           // Cumulative eta
        float prob;                 // Probability to continue
        bool foundIntersection;     // Whether there is an intersection point found

        /* Start the path tracing */
        for (uint32_t bounces = 0; ; ++bounces) {
            /* If there is no intersection found, break from the loop */
            foundIntersection = scene->rayIntersect(iterRay, its);
            if (!foundIntersection) break;

            /* Add the emitter light at path vertex */
            if (its.mesh->isEmitter()) {
                eQ = EmitterQueryRecord(-iterRay.d, its.shFrame.n);
                L += throughput * its.mesh->getEmitter()->eval(eQ);
            }
            
            /* We sample a new direction on this surface */
            BSDFQueryRecord bQ(its.shFrame.toLocal(-iterRay.d).normalized());
            bsdf = its.mesh->getBSDF();
            fr = bsdf->sample(bQ, sampler->next2D());
            if (fr.getLuminance() == 0.0f) break;

            /* Update the iteration ray given the sampling */
            Ray3f temp = Ray3f(its.p, its.shFrame.toWorld(bQ.wo));
            memcpy(&iterRay, &temp, sizeof(Ray3f));

            /* Compute the probability to continue */
            eta *= bQ.eta;
            prob = std::min(0.99f, throughput.maxCoeff() * eta * eta);
            /* Only start doing Russian Roulette after at least three bounces */
            prob = bounces < 4 ? 1 : prob;

            /* Use the Russian Roulette */
            if (sampler->next1D() < prob) {
                throughput *= fr / prob;
            } else {
                break;
            }
        }
        return L;
    }

    std::string toString() const {
        return "PathMatsIntegrator[]";
    }

};

NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats");
NORI_NAMESPACE_END