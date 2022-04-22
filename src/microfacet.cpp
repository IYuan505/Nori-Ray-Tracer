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

NORI_NAMESPACE_BEGIN

class Microfacet : public BSDF {
public:
    Microfacet(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);

        /* Albedo of the diffuse base material (a.k.a "kd") */
        m_kd = propList.getColor("kd", Color3f(0.5f));

        /* To ensure energy conservation, we must scale the 
           specular component by 1-kd. 

           While that is not a particularly realistic model of what 
           happens in reality, this will greatly simplify the 
           implementation. Please see the course staff if you're 
           interested in implementing a more realistic version 
           of this BRDF. */
        m_ks = 1 - m_kd.maxCoeff();
    }

    /// Evaluate the BRDF for the given pair of directions
    Color3f eval(const BSDFQueryRecord &bRec) const {
    	Color3f color_diffuse;           // color component from the diffuse
        Color3f color_dielectric;        // color component from the dielectric microfacet
        Vector3f wh;                     // middle vector between wi and wo
        float microfacet_dist;           // microfacet distributioin
        float fr;                        // fresnel coefficient
        float shadow_mask;               // shadowing and mask
        float gi, go;                    // shadowing and mask individually
        float bi, bo, ci, co;            // polynomial proximation of the shadowing and mask

        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        /* Compute the diffuse component. */
        color_diffuse = m_kd * INV_PI;

        /* Compute the dielectric microfacet component. 
           First, compute the microfacet distribution. */
        wh = (bRec.wi + bRec.wo).normalized();
        microfacet_dist = Warp::squareToBeckmannPdf(wh, m_alpha);

        /* Second, compute the fresnel ratio */
        fr = fresnel(wh.dot(bRec.wi), m_extIOR, m_intIOR);

        /* Third, compute the shadowing and mask */
        bi = 1 / (m_alpha * Frame::tanTheta(bRec.wi));
        ci = bRec.wi.dot(wh) / bRec.wi.z() > 0;
        gi = bi < 1.6 ? 
            ci * (3.535 * bi + 2.181 * bi * bi) / (1 + 2.276 * bi + 2.577 * bi * bi) : ci;
        bo = 1 / (m_alpha * Frame::tanTheta(bRec.wo));
        co = bRec.wo.dot(wh) / bRec.wo.z() > 0;
        go = bo < 1.6 ?
            co * (3.535 * bo + 2.181 * bo * bo) / (1 + 2.276 * bo + 2.577 * bo * bo) : co;
        shadow_mask = gi * go;


        /* Last, compute the normalization of dielectric microfacet component */
        color_dielectric = m_ks * microfacet_dist * fr * shadow_mask / 
            (4 * Frame::cosTheta(bRec.wi) * Frame::cosTheta(bRec.wo) * Frame::cosTheta(wh));

        return color_diffuse + color_dielectric;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
    	Vector3f wh;                     // middle vector between wi and wo
        float microfacet_dist;           // microfacet distributioin
        float jacobian;                  // Jacobian of the half direction mapping

        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        wh = (bRec.wi + bRec.wo).normalized();
        microfacet_dist = Warp::squareToBeckmannPdf(wh, m_alpha);
        jacobian = 1 / (4 * wh.dot(bRec.wo));

        return m_ks * microfacet_dist * jacobian 
            + (1 - m_ks) * INV_PI * Frame::cosTheta(bRec.wo);
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        Vector3f wh; // middle vector between wi and wo
        
        bRec.measure = ESolidAngle;
        bRec.eta = 1.0f;
    	if (_sample.x() > m_ks){
            /* Sample on diffuse reflection */
            bRec.wo = Warp::squareToCosineHemisphere(
                Point2f((_sample.x() - m_ks) / (1 - m_ks), _sample.y()));
        }
        else {
            /* Sample on specular reflection */
            wh = Warp::squareToBeckmann(
                Point2f(_sample.x() / m_ks, _sample.y()), m_alpha);
            bRec.wo = 2 * (bRec.wi.dot(wh)) * wh - bRec.wi;
        }

        if (bRec.measure != ESolidAngle
            || Frame::cosTheta(bRec.wi) <= 0
            || Frame::cosTheta(bRec.wo) <= 0)
            return 0.0f;

        return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return true;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f,\n"
            "  kd = %s,\n"
            "  ks = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR,
            m_kd.toString(),
            m_ks
        );
    }
private:
    float m_alpha;
    float m_intIOR, m_extIOR;
    float m_ks;
    Color3f m_kd;
};

NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END
