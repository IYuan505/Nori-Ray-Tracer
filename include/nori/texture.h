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
    Bitmap *readPNG2Bitmap(std::string file_path) const {
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
            }
        }
        return bitmap;
    }

    Normal3f perturbNormal(const Point2f &uv, const Normal3f shadingNormal) const {
        if (m_bump == nullptr)
            return shadingNormal;

        int y = uv.y() * m_bump->rows();
        int x = uv.x() * m_bump->cols();
        y = (int) clamp(y, 0, m_bump->rows()-1);
        x = (int) clamp(x, 0, m_bump->cols()-1);

        Vector3f perturbedNormal; 
        for (int i = 0; i < 3; ++i)
            perturbedNormal[i] = m_bump->coeffRef(y, x)[i];

        Frame f = Frame(shadingNormal);
        return f.toWorld(perturbedNormal);

    }

    /// Return the texture given the intersection uv
    virtual Color3f eval(const Point2f &uv) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/BSDF/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return ETexture; }

protected:
    Bitmap *m_bump = nullptr;
};

NORI_NAMESPACE_END
