#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>
#include <nori/emitter.h>
#include <nori/bsdf.h>

NORI_NAMESPACE_BEGIN

class EmitterWhittedIntegrator : public Integrator {
public:
    EmitterWhittedIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;           // Intersection information
        Intersection itsShadow;     // Intersection for shadow ray
        Ray3f shadowRay;            // A shadow ray from intersection point to the light
        Color3f emitterLight = {0}; // Light of emitter
        Color3f finalLight = {0};   // final rendered color
        Color3f fr;                 // Reflection or refraction factor
        const BSDF *bsdf;           // BSDF of the current mesh
        Emitter * emitter;          // The sampled emitter
        EmitterQueryRecord eQ;      // The emitter query record
        Vector3f pToLight;          // Vector from point to light
        float geo;                  // Geometric term
        int visibility;             // Visibility from p to light

        if (!scene->rayIntersect(ray, its)){
            return Color3f(0.0f);
        } else {
            bsdf = its.mesh->getBSDF();
            if (bsdf->isDiffuse()){
                /* Check the intersection point is a normal scene or emitter */
                if (its.mesh->isEmitter()) {
                    eQ = EmitterQueryRecord(-ray.d, its.shFrame.n);
                    finalLight = its.mesh->getEmitter()->eval(eQ);
                }
                
                emitter = scene->getEmitter();
                /* Sample points on the emitter */
                eQ = EmitterQueryRecord();
                eQ.itsP = its.p;
                emitter->sample(eQ, sampler);
                pToLight =  eQ.p - its.p;

                /* Compute the visibitlity */
                shadowRay.o = its.p;
                shadowRay.d = pToLight.normalized();
                /* Epsilon for shadow ray 1e-4 (important) */
                shadowRay.mint = 1e-4f;
                shadowRay.maxt = pToLight.norm() - 1e-4f;
                shadowRay.update();
                visibility = 1 - (int) scene->getAccel()->rayIntersect(shadowRay, itsShadow, true);
                /* Test the ray is on the correct side */
                visibility = visibility && pToLight.dot(its.shFrame.n) > 0 && pToLight.dot(eQ.n) < 0;
                if (visibility == 0)
                    return finalLight;

                /*  Compute the radiance */
                eQ.wo = - pToLight.normalized();
                emitterLight = emitter->eval(eQ);

                /* Compute the bsdf's fr */
                fr = bsdf->eval(
                    BSDFQueryRecord(its.shFrame.toLocal(pToLight), its.shFrame.toLocal(- ray.d), ESolidAngle));

                /* Compute the geometric term */
                geo = ((its.shFrame.n.normalized()).dot(pToLight.normalized()))
                    * ((eQ.n.normalized()).dot(- pToLight.normalized())) / (pToLight.dot(pToLight));
                geo = abs(geo);

                /* Final light, also considering emitter pdf */
                finalLight += Color3f(fr * emitterLight * geo / eQ.pdf);
                return finalLight;
            } else {
                /* Generate a reflection or refraction sample */
                BSDFQueryRecord bQ(its.shFrame.toLocal(- ray.d));
                fr = bsdf->sample(bQ, sampler->next2D());

                if (sampler->next1D() < 0.95) {
                    return fr * Li(scene, sampler, Ray3f(its.p, its.shFrame.toWorld(bQ.wo))) / 0.95;
                } else {
                    return 0.0f;
                }
            }
        }
    }

    std::string toString() const {
        return "EmitterWhittedIntegrator[]";
    }

};

NORI_REGISTER_CLASS(EmitterWhittedIntegrator, "emitterwhitted");
NORI_NAMESPACE_END