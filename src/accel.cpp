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

#include <nori/accel.h>
#include <nori/timer.h>
#include <tbb/parallel_for.h>
#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

// uint32_t tot_size = 0;
// uint32_t leaf_node = 0;
// uint32_t internal_node = 0;

void Accel::addMesh(Mesh *mesh) {
    if (m_mesh)
        throw NoriException("Accel: only a single mesh is supported!");
    m_mesh = mesh;
    m_bbox = m_mesh->getBoundingBox();
}

void Accel::build() {
    Timer timer;

    /* Set up the root OctreeNode */
    memcpy(&root.bbox, &m_bbox, sizeof(m_bbox));

    /* Initialize parameters to build the Octree */
    uint32_t size = m_mesh->getTriangleCount();
    uint32_t *triangles = (uint32_t *) malloc(sizeof(uint32_t) * size);
    for (uint32_t i = 0; i < size; ++i) {
        triangles[i] = i;
    }

    /* Build the Octree, with depth = 8 */
    // buildOctree(&root, triangles, size, 8);

    /* Build the Octree, with parallelism */
    buildOctree(&root, triangles, size, 1);

    tbb::parallel_for( tbb::blocked_range<int>(0, 8),
                       [&](tbb::blocked_range<int> r)
    {
        for (int i=r.begin(); i<r.end(); ++i)
        {
            buildOctree(root.children[i], root.children[i]->leafTriangles, 
                root.children[i]->leafSize, 7);
            root.children[i]->leafSize = (uint32_t) -1;
        }
    });

    cout << "Build Tree. (took " << timer.elapsedString() << ")" << endl;
    // printf("Number of leaf nodes: %d\nNumber of internal nodes: %d\nAverage triangles per leaf node: %f\nTot triangles: %d\n",
    //    leaf_node, internal_node, (float) tot_size / leaf_node, tot_size );
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its, bool shadowRay) const {
    bool foundIntersection = false;  // Was an intersection found so far?
    uint32_t f = (uint32_t) -1;      // Triangle index of the closest intersection

    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its '.maxt' value)

    /* Brute force search through all triangles */
    // for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
    //     float u, v, t;
    //     if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
            /* An intersection was found! Can terminate
               immediately if this is a shadow ray query */
    //         if (shadowRay)
    //             return true;
    //         ray.maxt = its.t = t;
    //         its.uv = Point2f(u, v);
    //         its.mesh = m_mesh;
    //         f = idx;
    //         foundIntersection = true;
    //     }
    // }

    /* Search using Octree data structure */
    if (root.bbox.rayIntersect(ray))
        searchOctree(&root, ray, its, f, shadowRay);
    foundIntersection = (f != (uint32_t) -1);

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1-its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh   = its.mesh;
        const MatrixXf &V  = mesh->getVertexPositions();
        const MatrixXf &N  = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F  = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) +
                bary.y() * UV.col(idx1) +
                bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1-p0).cross(p2-p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame = Frame(
                (bary.x() * N.col(idx0) +
                 bary.y() * N.col(idx1) +
                 bary.z() * N.col(idx2)).normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

void Accel::buildOctree(OctreeNode *currentNode, uint32_t *currentTriangles, 
    uint32_t size, uint32_t remainingDepth) {

    /* Local variables are only used when currentNode becomes an internal node */
    uint32_t *childrenTriangles[8]; // clean up is left for children
    uint32_t childrenTrianglesSize[8];
    uint32_t childrenTrianglesCapacity[8];
    uint32_t initSize;

    if (size < 10 || remainingDepth == 0 ) {
        /* If the size of the current triangles is small, make a leaf node or the 
         * remaining depth is 0, make a leaf node */
        currentNode->leafSize = size;
        currentNode->leafTriangles = currentTriangles;
        // tot_size += size;
        // leaf_node++;
    } else {
        // internal_node++;
        /* If neither conditions is satisfied, make 8 children OctreeNodes */
        currentNode->leafSize  = (uint32_t) -1; // indicates it is an internal node
        Point3f center = currentNode->bbox.getCenter();
        initSize = size / 8;
        for (uint32_t i = 0; i < 8; ++i) {
            /* Set up children OctreeNode */
            currentNode->children[i] = (OctreeNode *) malloc(sizeof(OctreeNode));
            currentNode->children[i]->bbox.reset();
            currentNode->children[i]->bbox.expandBy(currentNode->bbox.getCorner(i));
            currentNode->children[i]->bbox.expandBy(center);

            /* Prepare the temporary data structure that is used for \ref buildOctree */
            childrenTriangles[i] = (uint32_t *) malloc(sizeof(uint32_t) * initSize);
            childrenTrianglesSize[i] = 0;
            childrenTrianglesCapacity[i] = initSize;
        }

        /* Assign each triangle into different children OctreeNode */
        for (uint32_t triangleIdx = 0; triangleIdx < size; ++triangleIdx) {
            for (uint32_t childrenIdx = 0; childrenIdx < 8; ++childrenIdx) {
                /* If the triangle overlaps with the children Octree, put its index into 
                 * childrenTriangles[j] */
                if (currentNode->children[childrenIdx]->bbox.overlaps(
                    m_mesh->getBoundingBox(currentTriangles[triangleIdx]), false)) {
                    /* Handle the size of the array exceeding its capacity */
                    if (childrenTrianglesSize[childrenIdx] == 
                        childrenTrianglesCapacity[childrenIdx]) {
                        childrenTrianglesCapacity[childrenIdx] <<= 1;
                        childrenTriangles[childrenIdx] = 
                            (uint32_t *) realloc(childrenTriangles[childrenIdx], 
                            sizeof(uint32_t) * childrenTrianglesCapacity[childrenIdx]);
                    }
                    childrenTriangles[childrenIdx][childrenTrianglesSize[childrenIdx]] 
                        = currentTriangles[triangleIdx];
                    childrenTrianglesSize[childrenIdx]++;
                }
            }
        }

        free(currentTriangles);
        for (uint32_t i = 0; i < 8; ++i) {
            if (childrenTrianglesSize[i] > 0) {
                buildOctree(currentNode->children[i], childrenTriangles[i], 
                    childrenTrianglesSize[i], remainingDepth - 1);
            }
            else {
                /* Directly make a leaf, without recursion */
                currentNode->children[i]->leafSize = 0;
                // leaf_node++;
            }
        }
    }
}

void Accel::searchOctree(const OctreeNode *currentNode, Ray3f &ray, Intersection &its, 
    uint32_t &f, bool shadowRay) const { 

    float u, v, t, dummy;
    std::pair <int, float> childrenDistance[8]; /* First column: children's index
                                                 * Second column: children's nearT */

    if (currentNode->leafSize == (uint32_t) -1) {
        /* Unordered version of the ray traversal */
        // for (uint32_t i = 0; i < 8; ++i) {
        //     if (currentNode->children[i]->leafSize != 0
        //         && currentNode->children[i]->bbox.rayIntersect(ray)) {
        //         searchOctree(currentNode->children[i], ray, its, f, shadowRay);
        //         if (shadowRay && (f != (uint32_t) -1))
        //             return;
        //     }
        // }
        /* Ordered version of the ray traversal */
        for (uint32_t i = 0; i < 8; ++i) {
            childrenDistance[i].first = i;
            /* non empty children */
            if (currentNode->children[i]->leafSize != 0) {
                if (!currentNode->children[i]->bbox.rayIntersect(
                    ray, childrenDistance[i].second, dummy)) {
                    childrenDistance[i].second = std::numeric_limits<float>::infinity();
                }
            } else {
                childrenDistance[i].second = std::numeric_limits<float>::infinity();
            }
        }

        /* Lambda function that is used for sorting the children distance, based on nearT */
        std::sort(childrenDistance, childrenDistance + 8, 
            [](std::pair<int, float>& a, std::pair<int, float>& b){ return a.second < b.second;});

        for (uint32_t i = 0; i < 8; ++i) {
            if (currentNode->children[childrenDistance[i].first]->leafSize != 0 
                && childrenDistance[i].second != std::numeric_limits<float>::infinity()) {
                searchOctree(currentNode->children[childrenDistance[i].first], ray, its, f, shadowRay);
                /* If found f, and current t is smaller than next nearT, or it is a 
                 * shadaw ray, return */
                if ((f != (uint32_t) -1) && 
                    ((i < 7 && its.t <= childrenDistance[i+1].second) 
                        || shadowRay))
                    return;
            }
        }
    } else {
        for (uint32_t i = 0; i < currentNode->leafSize; ++i) {
            if (m_mesh->rayIntersect(currentNode->leafTriangles[i], ray, u, v, t)) {
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = m_mesh;
                f = currentNode->leafTriangles[i];
                if (shadowRay) return;
            }
        }
    }
}

NORI_NAMESPACE_END

