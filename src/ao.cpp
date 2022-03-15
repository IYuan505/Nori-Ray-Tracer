#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AoIntegrator : public Integrator {
public:
    AoIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;     // Intersection information
        Ray3f shadowRay;      // A shadow ray from intersection point to a sampled direction
        Vector3f sphereSample;// The sampled direction of the sphere

        if (!scene->rayIntersect(ray, its)){
            return Color3f(0.0f);
        } else {
            /* Sample points on the sphere with cosine distribution and convert the points
             * from the local coordinate system to the world coordinate system. The world
             * coordinate system is what we use for the ray direction */
            sphereSample = its.shFrame.toWorld(
                Warp::squareToCosineHemisphere(sampler->next2D()));
            shadowRay.o = its.p;
            shadowRay.d = sphereSample;
            shadowRay.update();
            return Color3f(1 - (int) scene->getAccel()->rayIntersect(shadowRay, its, true));
        }
    }

    std::string toString() const {
        return "AoIntegrator[]";
    }

};

NORI_REGISTER_CLASS(AoIntegrator, "ao");
NORI_NAMESPACE_END