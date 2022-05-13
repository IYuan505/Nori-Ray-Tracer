#include <nori/texture.h>
#include <nori/color.h>
#include <nori/bitmap.h>
#include <unistd.h>
#include <lodepng.h>

NORI_NAMESPACE_BEGIN

class ImageTexture: public Texture {
public:
    ImageTexture(const PropertyList &propList) {
        std::string texture_path = propList.getString("texture_path", "");
        if (access(texture_path.c_str(), F_OK ) != -1) {
            unsigned char *rgb;
            unsigned w, h;
            unsigned int error = lodepng_decode24_file(&rgb, &w, &h, texture_path.c_str());
            if (error != 0)
                throw NoriException("ImageTexture: unable to read PNG file!");
            
            m_bitmap = new Bitmap(Vector2i(w, h));
            for (unsigned int y = 0; y < h; ++y) {
                for (unsigned int x = 0; x < w; ++x, rgb += 3) {
                    m_bitmap->coeffRef(h - y - 1, x)[0] = rgb[0] / 255.f;
                    m_bitmap->coeffRef(h - y - 1, x)[1] = rgb[1] / 255.f;
                    m_bitmap->coeffRef(h - y - 1, x)[2] = rgb[2] / 255.f;
                }
            }
        } else {
            throw NoriException(
                "ImageTexture: image file does not exist!");
        }
    }

    Color3f eval(const Point2f &uv) const {
        int y = uv.y() * m_bitmap->rows();
        int x = uv.x() * m_bitmap->cols();
    	return m_bitmap->coeffRef(y, x);
    }

    /// Return a human-readable summary
    std::string toString() const {
        return "ImageTexture[]";
    }

private:
	Bitmap *m_bitmap = nullptr;
};

NORI_REGISTER_CLASS(ImageTexture, "imagetexture");
NORI_NAMESPACE_END
