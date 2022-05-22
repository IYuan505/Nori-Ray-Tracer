#include <nori/texture.h>
#include <nori/color.h>
#include <unistd.h>

NORI_NAMESPACE_BEGIN

class ConstantTexture: public Texture {
public:
    ConstantTexture(const PropertyList &propList) {
        m_value = propList.getColor("value", Color3f(1.0f));
        std::string bump_path = propList.getString("bump_path", "");
        /* Read bump map if there is any */
        if (access(bump_path.c_str(), F_OK ) != -1)
            m_bump = readPNG2Bitmap(bump_path);
    }

    Color3f eval(const Point2f &uv) const {
    	return m_value;
    }

    /// Return a human-readable summary
    std::string toString() const {
        return tfm::format(
            "ConstantTexture[\n"
            "  value = %s\n"
            "]", m_value.toString());
    }

private:
	Color3f m_value;
};

NORI_REGISTER_CLASS(ConstantTexture, "constanttexture");
NORI_NAMESPACE_END
