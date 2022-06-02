#include <nori/texture.h>
#include <nori/color.h>
#include <nori/bitmap.h>
#include <nori/frame.h>
#include <unistd.h>

NORI_NAMESPACE_BEGIN

class ImageTexture: public Texture {
public:
    ImageTexture(const PropertyList &propList) {
        std::string texture_path = propList.getString("texture_path", "");
        std::string bump_path = propList.getString("bump_path", "");

        /* Read image texture */
        if (access(texture_path.c_str(), F_OK ) != -1)
            m_bitmap = readPNG2Bitmap(texture_path, false);
        else 
            throw NoriException(
                "ImageTexture: image file does not exist!");
    }

    Color3f eval(const Point2f &uv) const {
        int y = uv.y() * m_bitmap->rows();
        int x = uv.x() * m_bitmap->cols();
        y = (int) clamp(y, 0, m_bitmap->rows()-1);
        x = (int) clamp(x, 0, m_bitmap->cols()-1);
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
