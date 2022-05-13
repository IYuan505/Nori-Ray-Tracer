#pragma once

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Texture for objects
 *
 * The texture class provides texture to objects
 */
class Texture: public NoriObject {
public:
    /// Return the texture given the intersection point
    virtual Color3f eval(const Point2f &uv) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return ETexture; }
};

NORI_NAMESPACE_END
