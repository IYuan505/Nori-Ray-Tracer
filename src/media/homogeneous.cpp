#include <nori/medium.h>

NORI_NAMESPACE_BEGIN

class HomogeneousMedium: public Medium {
public:
	HomogeneousMedium(const PropertyList &propList) {
		m_sigma_a = propList.getColor("sigma_a", Color3f(0.2f));
		m_sigma_s = propList.getColor("sigma_s", Color3f(0.3f));
		m_sigma_t = m_sigma_a + m_sigma_s;
		m_phase = new HenyeyGreenstein(propList.getFloat("g", 0.0f));
	}

	Color3f tr(const Ray3f ray, Sampler *sampler) const {
		float d = std::min(ray.maxt * ray.d.norm(), std::numeric_limits<float>::max());
	    return Color3f(std::exp(-m_sigma_t[0] * d), 
	    	std::exp(-m_sigma_t[1] * d), std::exp(-m_sigma_t[2] * d));
	}

	Color3f sample(const Ray3f ray, Sampler *sampler, Intersection &mediumIts) const {
	    // Sample a channel and distance along the ray
	    int channel = std::min((int)(sampler->next1D() * 3), 2);
	    float dist = -std::log(1 - sampler->next1D()) / m_sigma_t[channel];
	    float t = std::min(dist / ray.d.norm(), ray.maxt);
	    bool sampledMedium = t < ray.maxt;
	    
	    mediumIts.p = ray.o + ray.d * t;
	    mediumIts.t = t;

	    // Compute the transmittance and sampling density
	    float d = std::min(t, std::numeric_limits<float>::max()) * ray.d.norm();
	    Color3f tr = Color3f(std::exp(-m_sigma_t[0] * d), 
	    	std::exp(-m_sigma_t[1] * d), std::exp(-m_sigma_t[2] * d));

	    // Return weighting factor for scattering from homogeneous medium
	    Color3f density = sampledMedium ? (m_sigma_t * tr) : tr;
	    float pdf = (density[0] + density[1] + density[2]) / 3;
	    if (pdf == 0) return 0.0f;

	    tr = Color3f(tr[0] / pdf, tr[1] / pdf, tr[2] / pdf);
	    return sampledMedium ? (tr * m_sigma_s) : tr;
	}

    float sample_p(const Vector3f &wi, Vector3f *wo, const Point2f &u) const {
    	return m_phase->sample_p(wi, wo, u);
	}

	float p(const Vector3f &wi, const Vector3f &wo) const {
		return m_phase->p(wi, wo);
	}

	Color3f albedo() const {
		return m_sigma_s /  m_sigma_t;
	}

	Color3f getS() const {
		return m_sigma_s;
	}

	std::string toString() const {
        return "Homogeneous medium\n";
    }

private:
	Color3f m_sigma_a, m_sigma_s, m_sigma_t;
	const HenyeyGreenstein *m_phase;
};

NORI_REGISTER_CLASS(HomogeneousMedium, "homo");
NORI_NAMESPACE_END