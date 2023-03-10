#include <nori/bsdf.h>
#include <nori/medium.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Null BRDF model
 */
class NullBSDF : public BSDF {
public:
    NullBSDF(const PropertyList &propList) {}

    /// Evaluate the BRDF model
    Color3f eval(const BSDFQueryRecord &bRec) const { return 1.0f; }

    /// Compute the density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const { return 0.0f; }

    /// Draw a a sample from the BRDF model
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const { 
        bRec.measure = ESolidAngle;
        bRec.wo = -bRec.wi;
        bRec.eta = 1.0f;
        return 1.0f;
    }

    bool isNull() const {
        return true;
    }

    /// Return a human-readable summary
    std::string toString() const {
        return "Null BSDF";
    }

    EClassType getClassType() const { return EBSDF; }

};

NORI_REGISTER_CLASS(NullBSDF, "null");
NORI_NAMESPACE_END
