#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AreaEmitter : public Emitter {
public:
    AreaEmitter(const PropertyList &props) {
        radiance = props.getColor("radiance");
    }

    Color3f eval(const EmitterQueryRecord &eRec) const {
        /* Only return radiance when it is above the surface */
        float cosTheta = eRec.wo.dot(eRec.n);
        return cosTheta > 0 ? radiance : Color3f(0.0f);
    }

    void sample(EmitterQueryRecord &eRec, Sampler *sampler) const {
        /* Sampled points on the mesh */
        AreaSample *areaSample;
        areaSample = eRec.mesh->getAreaLightSample(sampler->next2D());
        eRec.p = areaSample->p;
        eRec.n = areaSample->n;
        eRec.pdf = areaSample->pdf;
    }

    float pdf(const EmitterQueryRecord &eRec) const {
        return eRec.pdf;
    }

    std::string toString() const {
        return tfm::format(
            "AreaEmitter[\n"
            "  radiance = %s\n]", 
            radiance.toString());
    }

private:
    Color3f radiance; // radience of the emitter
};

NORI_REGISTER_CLASS(AreaEmitter, "area");
NORI_NAMESPACE_END