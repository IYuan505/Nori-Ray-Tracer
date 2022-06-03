/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN

class RoughDielectric : public BSDF {
public:
    RoughDielectric(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }


    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
        float color_reflection;          // reflection color
        float color_refraction;          // refraction color
        Vector3f wh;                     // middle vector between wi and wo
        float microfacet_dist;           // microfacet distributioin
        float fr;                        // fresnel coefficient
        float shadow_mask;               // shadowing and mask
        float gi, go;                    // shadowing and mask individually
        float bi, bo, ci, co;            // polynomial proximation of the shadowing and mask
        float etaI, etaT, eta;           // eta for incident and transmission light, and ratio between two
        bool reflection;                 // whether it is reflection or refraction

        if (bRec.measure != ESolidAngle)
            return 0.0f;

        reflection = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) > 0;

        etaI = m_extIOR, etaT = m_intIOR;
        if (Frame::cosTheta(bRec.wi) < 0.0f)
            std::swap(etaI, etaT);
        eta = etaI / etaT;

        /* First, compute the microfacet distribution. */
        wh = reflection ? (bRec.wi + bRec.wo).normalized() : (eta * bRec.wi + bRec.wo).normalized();
        wh = wh.z() > 0 ? wh : -wh;
        microfacet_dist = Warp::squareToBeckmannPdf(wh, m_alpha);
        if (microfacet_dist == 0) return 0.0f;

        /* Second, compute the fresnel ratio */
        fr = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);

        /* Third, compute the shadowing and mask */
        bi = 1 / (m_alpha * Frame::tanTheta(bRec.wi));
        /* Check for NaN, tan(theta) = infinity */
        bi = std::isnan(bi) ? std::numeric_limits<float>::infinity() : bi;
        ci = bRec.wi.dot(wh) / bRec.wi.z() > 0;
        gi = bi < 1.6 ? 
            ci * (3.535 * bi + 2.181 * bi * bi) / (1 + 2.276 * bi + 2.577 * bi * bi) : ci;
        bo = 1 / (m_alpha * Frame::tanTheta(bRec.wo));
        bo = std::isnan(bo) ? std::numeric_limits<float>::infinity() : bo;
        co = bRec.wo.dot(wh) / bRec.wo.z() > 0;
        go = bo < 1.6 ?
            co * (3.535 * bo + 2.181 * bo * bo) / (1 + 2.276 * bo + 2.577 * bo * bo) : co;
        shadow_mask = gi * go;

        if (reflection) {
            /* Reflection */
            color_reflection = microfacet_dist * fr * shadow_mask / 
                (4 * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) * Frame::cosTheta(wh));
            return std::abs(color_reflection);
        } else {
            /* Refraction */
            color_refraction = std::abs(bRec.wi.dot(wh)) * std::abs(bRec.wo.dot(wh)) * microfacet_dist * (1 - fr) * shadow_mask /
                (- Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) * Frame::cosTheta(wh) * (eta * bRec.wi.dot(wh) + bRec.wo.dot(wh)) * (eta * bRec.wi.dot(wh) + bRec.wo.dot(wh)));
            return std::abs(color_refraction);
        }
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        Vector3f wh;                     // middle vector between wi and wo
        float microfacet_dist;           // microfacet distributioin
        float jacobian;                  // Jacobian of the half direction mapping
        float fr;                        // fresnel coefficient
        float etaI, etaT, eta;           // eta for incident and transmission light, and ratio between two
        float prob;                      // probability of sampling reflect/refract
        bool reflection;                 // whether it is reflection or refraction

        if (bRec.measure != ESolidAngle)
            return 0.0f;

        reflection = Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) > 0;

        etaI = m_extIOR, etaT = m_intIOR;
        if (Frame::cosTheta(bRec.wi) < 0.0f)
            std::swap(etaI, etaT);
        eta = etaI / etaT;

        wh = reflection ? (bRec.wi + bRec.wo).normalized() : - (eta * bRec.wi + bRec.wo).normalized();
        wh = wh.z() > 0 ? wh : -wh;
        fr = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
        prob = reflection ? fr : 1 - fr;
        jacobian = reflection ? std::abs(1 / (4 * wh.dot(bRec.wo))) : 
            std::abs(bRec.wo.dot(wh)) / ((eta * bRec.wi.dot(wh) + bRec.wo.dot(wh)) * (eta * bRec.wi.dot(wh) + bRec.wo.dot(wh)));
        microfacet_dist = Warp::squareToBeckmannPdf(wh, m_alpha);

        return prob * microfacet_dist * jacobian;
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        Vector3f wh;                     // middle vector between wi and wo
        Vector3f localWi;                // wi in local coordinate system
        Frame f;                         // frame for refraction computation
        float fr;                        // fresnel coefficient
        float etaI, etaT, eta;           // eta for incident and transmission light, and ratio between two
        float cosTheta, sinTheta;        // cos and sin theta of the incident light
        float sinThetaRefra;             // sin theta of the refraction light
        float cosThetaRefra;             // cos theta of the refraction light
        float zRefra;                    // z axis of the refracted light
        bool reflection;                 // sample reflection or refraction        

        bRec.measure = ESolidAngle;
        /* Reuse sample.x() */
        wh = Warp::squareToBeckmann(Point2f(_sample.x(), _sample.y()), m_alpha);
        fr = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);
        etaI = m_extIOR, etaT = m_intIOR;
        if (Frame::cosTheta(bRec.wi) < 0.0f)
            std::swap(etaI, etaT);
        eta = etaI / etaT;
        reflection = static_cast <float> (std::rand()) / static_cast <float> (RAND_MAX) <= fr;
        
        if (reflection){
            /* Sample on reflection */
            bRec.eta = 1.0f;
            bRec.wo = 2 * (bRec.wi.dot(wh)) * wh - bRec.wi;
            /* Make normal vector the same direction as wi */
            wh = Frame::cosTheta(bRec.wi) > 0 ? wh : -wh;
            if (bRec.measure != ESolidAngle
                || Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) <= 0
                || bRec.wi.dot(wh) <= 0)
                return 0.0f;
        } else {
            /* Sample on refraction */
            bRec.eta = bRec.wi.z() > 0.0f ? 
                m_intIOR / m_extIOR : m_extIOR / m_intIOR;

            /* Make normal vector the same direction as wi */
            wh = Frame::cosTheta(bRec.wi) > 0 ? wh : -wh;
            f = Frame(wh);
            localWi = f.toLocal(bRec.wi);
            cosTheta = std::min(1.0f, Frame::cosTheta(localWi));
            sinTheta = sqrt(1 - cosTheta * cosTheta);
            sinThetaRefra = sinTheta * eta;
            cosThetaRefra = sqrt(1 - sinThetaRefra * sinThetaRefra);
            zRefra = cosTheta > 0 ? - cosThetaRefra : cosThetaRefra;
            bRec.wo = f.toWorld(Vector3f(
                -localWi.x() * eta, -localWi.y() * eta, zRefra)).normalized();

            if (bRec.measure != ESolidAngle
                || Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) >= 0
                || bRec.wi.dot(wh) <= 0)
                return 0.0f;
        }

        return eval(bRec) * std::abs(Frame::cosTheta(bRec.wo)) / pdf(bRec);
    }


    std::string toString() const {
        return tfm::format(
            "RoughDielectric[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(RoughDielectric, "rough_dielectric");
NORI_NAMESPACE_END
