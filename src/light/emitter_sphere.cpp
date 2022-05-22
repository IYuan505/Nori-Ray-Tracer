#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class SphereEmitter : public Emitter {
public:
    SphereEmitter(const PropertyList &props) {
        radiance = props.getColor("radiance");
        position = props.getPoint("position");
        radius = props.getFloat("radius");
    }

    Color3f eval(const EmitterQueryRecord &eRec) const {
        /* Only return radiance when it is above the surface */
        float cosTheta = eRec.wo.dot(eRec.n);
        return cosTheta > 0 ? radiance : Color3f(0.0f);
    }

    void sample(EmitterQueryRecord &eRec, Sampler *sampler) const {
        /* Cosine hemisphere sampling */
        Point3f sample = Warp::squareToCosineHemisphere(sampler->next2D());
        eRec.pdf = Warp::squareToCosineHemispherePdf(sample) / (radius * radius);

        /* Rotate the hemishepre to the corret direction */
        Vector3f lightToP = (eRec.itsP - position).normalized();
        Vector3f x = Vector3f(lightToP.z(), lightToP.z(), - lightToP.x() - lightToP.y()).normalized();
        Vector3f y = Vector3f(lightToP.y()*x.z() - lightToP.z()*x.y(),
                            lightToP.z()*x.x() - lightToP.x()*x.z(),
                            lightToP.x()*x.y() - lightToP.y()*x.x()).normalized();
        
        Frame localFrame = Frame(x, y, lightToP);
        sample = localFrame.toWorld(sample);
        sample = sample * radius;
        
        eRec.p = Point3f(sample.x()+position.x(), sample.y()+position.y(), sample.z()+position.z());
        eRec.n = sample.normalized();
    }

    float pdf(const EmitterQueryRecord &eRec) const {
        return eRec.pdf;
    }

    std::string toString() const {
        return tfm::format(
            "SphereEmitter[\n"
            "  radiance = %s\n]", 
            radiance.toString());
    }

private:
    Color3f radiance;  // radience of the emitter
    Point3f position;  // center position of the sphere
    float radius;      // radius of the sphere

};

NORI_REGISTER_CLASS(SphereEmitter, "sphereEmitter");
NORI_NAMESPACE_END