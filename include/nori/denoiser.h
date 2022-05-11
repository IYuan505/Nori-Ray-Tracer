#pragma once

#include <nori/object.h>
#include <nori/block.h>
#include <nori/bitmap.h>
#include <nori/sampler.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Denoise the rendered image
 *
 * The denoise class provides denoising function for the bitmap
 */
class Denoiser: public NoriObject {
public:
    /// Denoise the bitmap
    virtual Bitmap *denoise(Bitmap *noisy, Sampler *sampler) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return EDenoiser; }
};

NORI_NAMESPACE_END
