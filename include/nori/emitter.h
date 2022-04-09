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

#pragma once

#include <nori/object.h>
#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Convenience data structure used to pass multiple
 * parameters to the evaluation and sampling routines in \ref Emitter
 */
struct EmitterQueryRecord {
    /// Outgoing direction (in the world frame)
    Vector3f wo;
    /// Position on the emitter
    Point3f p;
    /// Surface normal of the emitter
    Normal3f n;
    /// PDF of the sampling
    float pdf;
    /// Parent mesh of the emitter
    Mesh* mesh;
    /// Distance from its to emitter point
    float d;

    /// Improved version for sampling, also include intersection point info
    Point3f itsP;

    /// Empty init
    EmitterQueryRecord() {}

    /// Create a new record for sampling the Emitter
    EmitterQueryRecord(Mesh* mesh, Point3f itsP)
        : mesh(mesh), itsP(itsP) {}

    /// Create a new record for eval the Emitter
    EmitterQueryRecord(Vector3f wo, Normal3f n)
        : wo(wo), n(n) {}
};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
public:

    /**
     * \brief Evaluate the Emitter's radiance
     *
     * \param eRec
     *     A record with detailed information on the emitter query
     * 
     * \return
     *     The radiance value, evaluated for each color channel
     */
    virtual Color3f eval(const EmitterQueryRecord &eRec) const = 0;

    /**
     * \brief Evaluate the Emitter's radiance
     *
     * \param eRec
     *     A record with detailed information on the emitter query
     * \param sample  
     *     A uniformly distributed sample on \f$[0,1]^2\f$
     * 
     */
    virtual void sample(EmitterQueryRecord &eRec, Sampler *sampler) const = 0;

    /**
     * \brief Compute the probability of sampling \c eRec.p
     *
     * This method provides access to the probability density that
     * is realized by the \ref sample() method.
     *
     * \param eRec
     *     A record with detailed information on the emitter query
     *
     * \return
     *     A probability/density value expressed with respect
     *     to the specified measure
     */
    virtual float pdf(const EmitterQueryRecord &eRec) const = 0;
    
    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.) 
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }
};

NORI_NAMESPACE_END
