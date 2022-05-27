#pragma once

#include <nori/object.h>
#include <nori/frame.h>
#include <nori/mesh.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN


class PhaseFunction {
public:
    /* Returns the value of the phase function for the given pair of directions */
    virtual float p(const Vector3f &wo, const Vector3f &wi) const = 0;

    /* Get a incident direction for given outgoing direction, with additional random 2D u */
    virtual float sample_p(const Vector3f &wo, Vector3f *wi, const Point2f &u) const = 0;
};

/* Henyey–Greenstein phase function */
inline float PhaseHG(float cosTheta, float g) {
    float denom = 1 + g * g + 2 * g * cosTheta;
    return INV_FOURPI * (1 - g * g) / (denom * std::sqrt(denom));
}

class HenyeyGreenstein: public PhaseFunction {
public:
    HenyeyGreenstein(float g): g(g) {}

    float p(const Vector3f &wo, const Vector3f &wi) const {
        return PhaseHG(wo.dot(wi), g);
    }

    float sample_p(const Vector3f &wo, Vector3f *wi, const Point2f &u) const {
        float cosTheta;
        if (std::abs(g) < 1e-3)
            cosTheta = 1 - 2 * u[0];
        else {
            float sqrTerm = (1 - g * g) / (1 + g - 2 * g * u[0]);
            cosTheta = -(1 + g * g - sqrTerm * sqrTerm) / (2 * g);
        }

        // Compute direction _wi_ for Henyey--Greenstein sample
        float sinTheta = std::sqrt(std::max((float)0, 1 - cosTheta * cosTheta));
        float phi = 2 * M_PI * u[1];
        Frame f = Frame(wo);
        *wi = f.toWorld(Vector3f(sinTheta * std::cos(phi), sinTheta * std::sin(phi), cosTheta));
        return PhaseHG(cosTheta, g);
    }
private:
    const float g;
};

class Medium: public NoriObject {
public:
    /* Compute the beam transmittance */
    virtual Color3f tr(const Ray3f ray, Sampler *sampler) const = 0;
    /* Sample an intersection in the medium */
    virtual Color3f sample(const Ray3f ray, Sampler *sampler, Intersection &mediumIts) const = 0;
    /* Wrapper of the phase function */
    virtual float sample_p(const Vector3f &wo, Vector3f *wi, const Point2f &u) const = 0;
    virtual float p(const Vector3f &wo, const Vector3f &wi) const = 0;
    virtual Color3f albedo() const = 0;
    virtual Color3f getS() const = 0;
    EClassType getClassType() const { return EMedium; }
};

NORI_NAMESPACE_END
