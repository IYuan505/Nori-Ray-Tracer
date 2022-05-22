#include <nori/emitter.h>
#include <nori/scene.h>
#include <nori/sampler.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class RectangleEmitter : public Emitter {
public:
    RectangleEmitter(const PropertyList &props) {
        radiance = props.getColor("radiance");
        position = props.getPoint("position");
        width = props.getFloat("width");
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
        sample = Point2f((x + sampler->next1D()) * invSqrtVal * width - width/2,
                                (y + sampler->next1D()) * invSqrtVal * width - width/2);

        /* Only consider z facing square emitter */
        eRec.p = Point3f(sample.x(), position.y(), sample.y());
        eRec.n = n;
        eRec.pdf = 1 / (width * width);
    }

    float pdf(const EmitterQueryRecord &eRec) const {
        return eRec.pdf;
    }

    std::string toString() const {
        return tfm::format(
            "RectangleEmitter[\n"
            "  radiance = %s\n]", 
            radiance.toString());
    }

private:
    Color3f radiance;  // radience of the emitter
    Point3f position;  // center position of the rectangle
    Normal3f n;        // surface norm of the rectangle
    float width;       // width of the rectangle

};

NORI_REGISTER_CLASS(RectangleEmitter, "rectEmitter");
NORI_NAMESPACE_END