#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
public:
    SimpleIntegrator(const PropertyList &props) {
        lightPosition = props.getPoint("position");
        lightEnergy = props.getColor("energy");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Intersection its;     // Intersection information
        Ray3f shadowRay;      // A shadow ray from intersection point to the light
        Vector3f pToLight;    // A vector from intersection point to the light
        float cosTheta;       // The cosine of theta, where theta is the incidence angle
        float distance2;      // Distance square between light position and intersection point
        float factor;         // The factor of lightEnergy decay
        int visibility;       // If x and p are mutually visible

        if (!scene->rayIntersect(ray, its)){
            return Color3f(0.0f);
        } else {
            pToLight = lightPosition - its.p;
            distance2 = pToLight.dot(pToLight);
            cosTheta = pToLight.dot(its.shFrame.n) / (pToLight.norm() * its.shFrame.n.norm());
            shadowRay.o = its.p;
            shadowRay.d = pToLight;
            shadowRay.update();
            visibility = 1 - (int) scene->getAccel()->rayIntersect(shadowRay, its, true);
            factor = std::max(0.0f, cosTheta) * visibility / (distance2 * 4 * M_PI * M_PI);
            return Color3f(lightEnergy[0] * factor, lightEnergy[1] * factor, lightEnergy[2]*factor);
        }
    }

    std::string toString() const {
        return tfm::format(
            "SimpleIntegrator[\n"
            "  light position = %s,\n"
            "  energy = %s\n]", 
            lightPosition.toString(), lightEnergy.toString());
    }

private:
    Point3f lightPosition; // The position of the single point light
    Color3f lightEnergy;   // The energy color of the single point light
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END