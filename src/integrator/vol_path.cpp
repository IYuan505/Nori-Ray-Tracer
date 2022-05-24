#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>
#include <nori/vector.h>

NORI_NAMESPACE_BEGIN

class VolPathIntegrator : public Integrator {
public:
    VolPathIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;           // Intersection information
        Intersection mediumIts;     // Medium intersection
        Color3f L = {0};            // Final rendered light 
        Color3f throughput = {1};   // Path throughput
        Color3f fr;                 // Reflection or refraction factor (albedo)
        Ray3f iterRay(ray);         // Ray used for path tracing
        const Medium *currentMedium;// Current medium of the ray segment
        EmitterQueryRecord eQ;      // The emitter query record
        float eta = 1.0f;           // Cumulative eta
        float prob;                 // Probability to continue
        bool foundIntersection;     // Whether there is an intersection point found
        bool specularBounce = false;// Whether last bounce is a specular bounce

        /* Initialize current medium to the scene's global medium */
        currentMedium = scene->getMedium();

        /* Start the volumetric path tracing */
        for (uint32_t bounces = 0; ; ++bounces) {
            /* If there is no intersection found, break from the loop */
            foundIntersection = scene->rayIntersect(iterRay, its);
            iterRay.maxt = its.t;
            
            /* Sample the participating medium, if there is any */
            if (currentMedium != nullptr) throughput *= currentMedium->sample(iterRay, sampler, mediumIts);
            else mediumIts.t = std::numeric_limits<float>::infinity();

            if (throughput.getLuminance() == 0.0f) break;
            
            /* Intersect with a medium (true) or a surface (false) */
            if (mediumIts.t < its.t) {
                /* Handle scattering at point in medium for volumetric path tracer */
                L += throughput * scene->uniformlySampleLight(sampler, &mediumIts, &eQ, iterRay, currentMedium, false, true);

                Vector3f wo = -iterRay.d, wi;
                currentMedium->sample_p(wo, &wi, sampler->next2D());
                Ray3f temp = Ray3f(mediumIts.p, wi.normalized());
                memcpy(&iterRay, &temp, sizeof(Ray3f));
                specularBounce = false;
            } else {
                if (!foundIntersection) break;
                if ((bounces == 0 || specularBounce) && its.mesh->isEmitter()) {
                    eQ = EmitterQueryRecord(-iterRay.d, its.shFrame.n);
                    L += throughput * its.mesh->getEmitter()->eval(eQ);
                }

                /* Boundary of the medium object */
                if (its.mesh->getBSDF()->isNone()) {
                    Ray3f temp = Ray3f(its.p, iterRay.d);
                    memcpy(&iterRay, &temp, sizeof(Ray3f));
                    bounces--;

                    /* Enter into the medium */
                    if (currentMedium != its.mesh->getBSDF()->getIntMedium()) {
                        its.mesh->getBSDF()->setOutMedium(currentMedium);
                        currentMedium = its.mesh->getBSDF()->getIntMedium();
                    }
                    else {
                        currentMedium = its.mesh->getBSDF()->getOutMedium();
                    }
                    continue;
                }


                /* Explicit sampling direct emitters */
                L += throughput * scene->uniformlySampleLight(sampler, &its, &eQ, iterRay, currentMedium, true, true);

                /* Account for indirect light, we sample a new direction on this surface */
                BSDFQueryRecord bQ(its.shFrame.toLocal(-iterRay.d).normalized());
                fr = its.mesh->getBSDF()->sample(bQ, sampler->next2D());
                fr *= its.mesh->getTexture()->eval(its.uv);
                throughput *= fr;
                eta *= bQ.eta;
                if (fr.getLuminance() == 0.0f) break;

                /* Update the iteration ray given the sampling */
                Ray3f temp = Ray3f(its.p, its.shFrame.toWorld(bQ.wo).normalized());
                memcpy(&iterRay, &temp, sizeof(Ray3f));
                specularBounce = bQ.measure == EDiscrete;
            }

            /* Compute the probability to continue */
            prob = std::min(0.95f, throughput.maxCoeff() * eta * eta);
            /* Only start doing Russian Roulette after at least three bounces */
            prob = bounces < 4 ? 1 : prob;

            if (prob == 0) break;

            /* Use the Russian Roulette */
            if (sampler->next1D() < prob) {
                throughput /= prob;
            } else {
                break;
            }
        }
        return L;
    }

    std::string toString() const {
        return "VolPathIntegrator[]";
    }

};

NORI_REGISTER_CLASS(VolPathIntegrator, "vol_path");
NORI_NAMESPACE_END