#pragma once

#include <nori/object.h>
#include <nori/bitmap.h>
#include <nori/frame.h>
#include <lodepng.h>

NORI_NAMESPACE_BEGIN


/**
 * \brief Texture for objects
 *
 * The texture class provides texture to objects
 */
class Texture: public NoriObject {
public:
    Bitmap *readPNG2Bitmap(std::string file_path, bool raw) const {
        unsigned char *rgb;
        unsigned w, h;
        unsigned int error = lodepng_decode24_file(&rgb, &w, &h, file_path.c_str());
        if (error != 0)
            throw NoriException("ImageTexture: unable to read PNG file!");
        
        Bitmap *bitmap = new Bitmap(Vector2i(w, h));
        for (unsigned int y = 0; y < h; ++y) {
            for (unsigned int x = 0; x < w; ++x, rgb += 3) {
                bitmap->coeffRef(h - y - 1, x)[0] = clamp(rgb[0] / 255.f, 0.0f, 1.0f);
                bitmap->coeffRef(h - y - 1, x)[1] = clamp(rgb[1] / 255.f, 0.0f, 1.0f);
                bitmap->coeffRef(h - y - 1, x)[2] = clamp(rgb[2] / 255.f, 0.0f, 1.0f);
                if (!raw)
                    bitmap->coeffRef(h - y - 1, x) = bitmap->coeffRef(h - y - 1, x).toLinearRGB();
            }
        }
        return bitmap;
    }

    /// Return the texture given the intersection uv
    virtual Color3f eval(const Point2f &uv) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return ETexture; }

};

NORI_NAMESPACE_END
