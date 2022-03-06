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

#include <nori/mesh.h>

#define SELF_DEBUG true

NORI_NAMESPACE_BEGIN

/**
 * \brief Acceleration data structure used by \ref Accel
 */
struct OctreeNode {
    OctreeNode *children[8];         // The children OctreeNode if not leaf
    BoundingBox3f bbox;              // The bounding box of the OctreeNode
    uint32_t *leafTriangles;         // A leaf triangle array
    uint32_t leafSize;               /* size of the leaf triangle array, if it
                                      * is -1 (MAX_UINT32), then it is an 
                                      * internal node */
};

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {

public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its, bool shadowRay) const;

private:
    Mesh         *m_mesh = nullptr; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;           ///< Bounding box of the entire scene
    OctreeNode    root;             ///< Root of the Octree

    /** 
     * \brief Build the Octree data structure, called by \ref build()
     * 
     * \param currentNode
     *    A pointer points to the address of the current OctreeNode
     * 
     * \param currentTriangles
     *    A list of indexes of the triangles
     * 
     * \param size
     *    The size of the currentTriangles array
     * 
     * \param remainingDepth
     *    Remaining depth of the Octree to set up
     */
    void buildOctree(OctreeNode *currentNode, uint32_t *currentTriangles, uint32_t size, uint32_t remainingDepth);

    /**
     * \breif Find the closest hit of the ray to the triangles using Octree
     * 
     * \param currentNode
     *    A pointer points to the address of the current OctreeNode
     * 
     * \param ray
     *    The ray to search for hit
     * 
     * \param f
     *    The idx of the triangle with the closest hit
     * 
     * \param shadow
     *    true if this is a shadow ray query
     */
    void searchOctree(const OctreeNode *currentNode, Ray3f &ray, Intersection &its, 
        uint32_t &f, bool shadowRay) const;
};

NORI_NAMESPACE_END
