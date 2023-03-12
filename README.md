
# Nori Ray Tracer

Nori is an open-source, physically based renderer for generating realistic computer-generated images. It is primarily used in computer graphics and animation to create high-quality images of virtual environments and objects. Nori simulates the physics of light and its interaction with surfaces, which enables it to generate images with accurate lighting, shadows, reflections, and refractions. Nori is capable of producing photorealistic images with complex lighting and shading effects.

## 🐡 Feactured Results
A skull inside the homogenous participating media rendered using explicit emitter sampling, with image textures and depth-of-field enabled, denoised using BM3D denoiser.

![Rendered image with homogeneous participating media](results/project/project_final.png)

## File Structures
1. `results`: This folder includes some intermediate reports during the development process. It could be useful to understand how each individual feature is implemented.
2. `scenes`: This folder contains some scenes that could be rendered using the rendering engine. A typical setup would be: 1 XML file defining the scene (camera model, integrator, object positions, object materials...), and multiple object mesh files (defining objects themselves).
3. `src`: It is where the implementation actually resides. It is further divided into sub-directories by functionalities. 

## Features
### 🐛 Integrators
1. **Normal integrator**: returns the normal of the object in the view field.
2. **Point light integrator**: the light source is a point, and only direct illumination is considered.
3. **Ambient occlusion**: the light comes from a hemisphere, over the objects. It simulates the ambient lighting casting out onto your scene. Only direct illumination is considered.
4. **Simple integrator**: the difference between the simple integrator and the above ones is that the simple integrator can take area light as input. The light source can be from any position inside the scene.
5. **Whitted integrator**: whitted integrator fixed some issues with the previous integrators, by considering partial indirect lights. If the light hits the reflection or refraction surface, we will be continuing sampling, to consider the indirect illumination.
6. **Brute force path integrator**: the major difference between the path integrators and the above is that the indirect lights are all sampled with path integrators. The indirect lights are sampled using the surface properties. Lights can come from other surfaces of non-illuminating objects. A light path can be: emitter → surface1 → surface2 → ... → the point to be sampled.
7. **Path integrator with next event estimation**: compared to the brute force path integrator, it improves the efficiency of rendering by including explicit sampling light sources during the rendering time. This can be very useful if the light source is small, namely, the likelihood of directly sampling the surface to hit the light source is low. However, it can be a bad idea for non-diffuse materials, because the explicit emitter sampling is likely to produce no contribution.
8. **📍 Path integrator with multiple importance sampling**: MIS path integrator combines the advantages of the above two path integrators. It both sampling directions based on surface and directions from the emitters explicitly.
9. **📌 Volumetric path integrator with next event estimation**: this path integrator can render scenes with participating media inside, like fog, smoke, and so on. It uses next event estimation (explicit emitter sampling) to improve the rendering efficiency.  

### Materials
1. **Diffuse**: diffuse materials are the simplest materials that you could have. It reflects incoming lights in all directions at equal probability in the hemisphere where the surface normal points to.
2. **Dielectric**: dielectric materials resemble glass materials in real life. There can be both reflection and refraction when light hits this kind of object.
3. **Mirror**: as the name suggests, it implements a perfect mirror. 
4. **Microfacet**: this kind of material is most common in life. It has the roughness property, i.e., from smooth to rough.
5. **Rough dielectric**: it supports roughness on dielectric materials.
6. **Conductor**: normal conductors in life could be all metals, for example, gold, sodium, silver... The major difference between conductors and others is that conductors absorb and refract light.
7. **🕸 Null**: it is used for defined invisible media (fog, smoke, ...) boundaries.

### 🧤Denoiser
The noise of rendering could be reduced by increasing the sample count. However, to reduce the noise sigma to half, four times of sample count is required on average. The exponential increase makes it extremely hard to obtain a noise-free image for some scenes. Denoiser is a post-processing technic to remove noise from the image directly. Depending on the denoisng algorithm used, there could be some observable image quality degradation. This project implements a BM3D denoiser, which is the best-performant classic image denoiser (without using deep learning). 

### Cameras
1. **Perspective**: it is basically a pinhole camera model. It is an idealistic camera model without any distortion, or circle of confusion...
2. **Thin lens**: it is a more realistic camera model, supporting depth of field. With varying focal radii, the effects differ. With this camera model, you could realize the effects of depth-of-field.

### Light Sources
1. **Point light**: the light source is a point, emitting lights in all directions uniformly.
2. **Area light**: the light source is a 3D object, namely an object is emitting lights. If the 3D object is infinitely small, it is reduced to point lights. The area can be arbitrarily large as well. Note that the light is only cast on the side where the area norm points to. In other words, there is no light cast into the object's internals.

### Textures
The implementation supports image textures on objects. The mapping implemented is 2D _(u, v)_ mapping, which is the most widely used in real-world scenarios.
