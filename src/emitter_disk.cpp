#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class DiskEmitter : public Emitter {
public:
    DiskEmitter(const PropertyList &props) {
        radiance = props.getColor("radiance");
        position = props.getPoint("position");
        radius = props.getFloat("radius");
        n = props.getVector("norm");
    }

    Color3f eval(const EmitterQueryRecord &eRec) const {
        /* Only return radiance when it is above the surface */
        float cosTheta = eRec.wo.dot(eRec.n);
        return cosTheta > 0 ? radiance : Color3f(0.0f);
    }

    void sample(EmitterQueryRecord &eRec, Sampler *sampler) const {
        /* Stratified sampling */
        uint32_t sampleCnt = sampler->getSampleCount();
        int sqrtVal = (int) (std::sqrt((float) sampleCnt) + 0.5f);
        float invSqrtVal = 1.f / sqrtVal;
        int y = sampler->i / sqrtVal, x = sampler->i % sqrtVal;
        Point2f sample;
        sample = Point2f((x + sampler->next1D()) * invSqrtVal,
                                (y + sampler->next1D()) * invSqrtVal);
        sample = Warp::squareToUniformDisk(sample);
        sample = sample * radius;

        /* Only consider z facing square emitter */
        eRec.p = Point3f(sample.x(), position.y(), sample.y());
        eRec.n = n;
        eRec.pdf = 1 * INV_PI / (radius * radius);
    }

    float pdf(const EmitterQueryRecord &eRec) const {
        return eRec.pdf;
    }

    std::string toString() const {
        return tfm::format(
            "DiskEmitter[\n"
            "  radiance = %s\n]", 
            radiance.toString());
    }

private:
    Color3f radiance;  // radience of the emitter
    Point3f position;  // center position of the disk
    Normal3f n;        // surface norm of the disk
    float radius;           // radius of the disk

};

NORI_REGISTER_CLASS(DiskEmitter, "diskEmitter");
NORI_NAMESPACE_END