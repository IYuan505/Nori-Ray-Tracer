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

#include <nori/scene.h>
#include <nori/bitmap.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/camera.h>
#include <nori/emitter.h>
#include <nori/denoiser.h>

NORI_NAMESPACE_BEGIN

Scene::Scene(const PropertyList &) {
    m_accel = new Accel();
}

Scene::~Scene() {
    delete m_accel;
    delete m_sampler;
    delete m_camera;
    delete m_integrator;
    delete m_denoiser;
    delete m_medium;
}

void Scene::activate() {
    m_accel->build();

    if (!m_integrator)
        throw NoriException("No integrator was specified!");
    if (!m_camera)
        throw NoriException("No camera was specified!");
    
    if (!m_sampler) {
        /* Create a default (independent) sampler */
        m_sampler = static_cast<Sampler*>(
            NoriObjectFactory::createInstance("independent", PropertyList()));
    }

    cout << endl;
    cout << "Configuration: " << toString() << endl;
    cout << endl;
}

void Scene::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
        case EMesh: {
                Mesh *mesh = static_cast<Mesh *>(obj);
                m_accel->addMesh(mesh);
                m_meshes.push_back(mesh);
            }
            break;
        
        case EEmitter: {
                m_emitter = static_cast<Emitter *>(obj);
                /* TBD */
                // throw NoriException("Scene::addChild(): You need to implement this for emitters");
            }
            break;

        case ESampler:
            if (m_sampler)
                throw NoriException("There can only be one sampler per scene!");
            m_sampler = static_cast<Sampler *>(obj);
            break;

        case ECamera:
            if (m_camera)
                throw NoriException("There can only be one camera per scene!");
            m_camera = static_cast<Camera *>(obj);
            break;
        
        case EIntegrator:
            if (m_integrator)
                throw NoriException("There can only be one integrator per scene!");
            m_integrator = static_cast<Integrator *>(obj);
            break;

        case EDenoiser:
            if (m_denoiser)
                throw NoriException("There can only be one denoiser per scene!");
            m_denoiser = static_cast<Denoiser *>(obj);
            break;

        case EMedium:
            if (m_medium)
                throw NoriException("There can only be one global medium per scene!");
            m_medium = static_cast<Medium *>(obj);
            break;
            
        default:
            throw NoriException("Scene::addChild(<%s>) is not supported!",
                classTypeName(obj->getClassType()));
    }
}

Color3f Scene::uniformlySampleLight(Sampler *sampler, Intersection *its, 
    EmitterQueryRecord *eQ, const Ray3f &ray, const Medium *initMedium, bool surfaceIntersection) const {
    Intersection itsShadow;     // Intersection for shadow ray
    Ray3f shadowRay;            // A shadow ray from intersection point to the light
    const Medium *currentMedium;// Current medium of the ray
    Color3f emitterLight = {0}; // Light of emitter
    Color3f fr;                 // Reflection or refraction factor
    Color3f tr = {1.0f};        // Medium transmittance
    uint32_t chosenIdx;         // Chosen idx of emitter
    uint32_t emitterIdxSize = 0;// Number of emitters
    Mesh *mesh = NULL;          // The sampled emitter's mesh
    Emitter * emitter;          // The sampled emitter
    Vector3f pToLight;          // Vector from point to light
    float maxt;                 // Maxt of the shadow ray
    int visibility;             // Visibility from p to light

    currentMedium = initMedium;
    /* Explicit direct emitter sampling
       Select a sampling emitter uniformly */
    for (uint32_t i = 0; i < m_meshes.size(); ++i) {
        if (m_meshes[i]->isEmitter())
            emitterIdxSize++;
    }
    chosenIdx = sampler->next1D() * emitterIdxSize;
    for (uint32_t i = 0; i < m_meshes.size(); ++i) {
        if (m_meshes[i]->isEmitter()){
            if (chosenIdx == 0){
                mesh = m_meshes[i];
                break;
            }
            chosenIdx--;
        }
    }
    emitter = mesh->getEmitter();

    /* Sample points on the emitter */
    *eQ = EmitterQueryRecord(mesh, its->p);
    emitter->sample(*eQ, sampler);
    pToLight =  eQ->p - its->p;
    eQ->d = pToLight.norm();

    /* Compute the visibitlity, epsilon for shadow ray 1e-4 (important) */
    maxt = pToLight.norm() - 1e-4f;
    Ray3f temp = Ray3f(its->p, pToLight.normalized(), 1e-4f, maxt);
    memcpy(&shadowRay, &temp, sizeof(Ray3f));
    
    while (1) {
        visibility = 1 - (int) m_accel->rayIntersect(shadowRay, itsShadow, true);
        if (visibility == 0 && itsShadow.mesh->getBSDF()->isNull()) {
            shadowRay.maxt = itsShadow.t;
            if (currentMedium != nullptr)
                tr *= currentMedium->tr(shadowRay, sampler);

            if (itsShadow.shFrame.cosTheta(itsShadow.shFrame.toLocal(shadowRay.d)) > 0) {
                /* Coming from inside the medium */
                currentMedium = m_medium;
            } else {
                /* Comming outside from the medium */
                currentMedium = itsShadow.mesh->getBSDF()->getIntMedium();
            }

            maxt -= itsShadow.t;
            Ray3f temp = Ray3f(itsShadow.p, shadowRay.d, 1e-4f, maxt);
            memcpy(&shadowRay, &temp, sizeof(Ray3f));
        } else {
            if (visibility == 1 && currentMedium != nullptr) {
                /* No intersection is found, use the maxt of shadowRay directly */
                tr *= currentMedium->tr(shadowRay, sampler);
            }
            break;
        }
    }

    if (visibility == 0)
        return 0.0f;

    /*  Compute the radiance */
    eQ->wo = - pToLight.normalized();
    emitterLight = emitter->eval(*eQ);
    eQ->pdf *= pToLight.dot(pToLight) / abs(eQ->n.normalized().dot(-pToLight.normalized()));

    
    if (surfaceIntersection) {
        /* Compute the bsdf's fr */
        fr = its->mesh->getBSDF()->eval(
            BSDFQueryRecord(its->shFrame.toLocal(pToLight).normalized(), 
                its->shFrame.toLocal(-ray.d).normalized(), ESolidAngle));
        fr *= its->mesh->getTexture()->eval(its->uv);
        fr *= abs(its->shFrame.n.normalized().dot(pToLight.normalized()));
    } else {
        fr = initMedium->albedo() * initMedium->p(-ray.d.normalized(), pToLight.normalized());
    }

    return Color3f(fr * emitterLight * tr * emitterIdxSize / eQ->pdf);
}


std::string Scene::toString() const {
    std::string meshes;
    for (size_t i=0; i<m_meshes.size(); ++i) {
        meshes += std::string("  ") + indent(m_meshes[i]->toString(), 2);
        if (i + 1 < m_meshes.size())
            meshes += ",";
        meshes += "\n";
    }

    return tfm::format(
        "Scene[\n"
        "  integrator = %s,\n"
        "  sampler = %s\n"
        "  camera = %s,\n"
        "  meshes = {\n"
        "  %s  }\n"
        "]",
        indent(m_integrator->toString()),
        indent(m_sampler->toString()),
        indent(m_camera->toString()),
        indent(meshes, 2)
    );
}

NORI_REGISTER_CLASS(Scene, "scene");
NORI_NAMESPACE_END
