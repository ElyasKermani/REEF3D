/*--------------------------------------------------------------------
REEF3D
Copyright 2008-2025 Hans Bihs

This file is part of REEF3D.

REEF3D is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.
--------------------------------------------------------------------
Author: Elyas Larkermani
--------------------------------------------------------------------*/

#ifndef SIXDOF_BVH_H_
#define SIXDOF_BVH_H_

#include<Eigen/Dense>
#include<vector>
#include<algorithm>
#include<limits>
#include<cmath>

// Squared distance from p to the closest point on segment ab (Ericson)
inline double point_segment_sqdist(const Eigen::Vector3d& p,
                                   const Eigen::Vector3d& a,
                                   const Eigen::Vector3d& b)
{
    const Eigen::Vector3d ab = b - a;
    double t = ab.dot(p - a);
    const double denom = ab.squaredNorm();
    if(denom < 1.0e-30)
    {
        return (p - a).squaredNorm();
    }
    t = t / denom;
    if(t < 0.0) t = 0.0;
    else if(t > 1.0) t = 1.0;
    const Eigen::Vector3d q = a + t * ab;
    return (p - q).squaredNorm();
}

// Squared distance from p to triangle (a,b,c) — closest point in Voronoi/face (Ericson RTCD)
inline double point_triangle_sqdist(const Eigen::Vector3d& p,
                                    const Eigen::Vector3d& a,
                                    const Eigen::Vector3d& b,
                                    const Eigen::Vector3d& c)
{
    const Eigen::Vector3d ab = b - a;
    const Eigen::Vector3d ac = c - a;
    const Eigen::Vector3d ap = p - a;

    const double d1 = ab.dot(ap);
    const double d2 = ac.dot(ap);
    if(d1 <= 0.0 && d2 <= 0.0)
    {
        return ap.squaredNorm();
    }

    const Eigen::Vector3d bp = p - b;
    const double d3 = ab.dot(bp);
    const double d4 = ac.dot(bp);
    if(d3 >= 0.0 && d4 <= d3)
    {
        return bp.squaredNorm();
    }

    const Eigen::Vector3d cp = p - c;
    const double d5 = ab.dot(cp);
    const double d6 = ac.dot(cp);
    if(d6 >= 0.0 && d5 <= d6)
    {
        return cp.squaredNorm();
    }

    const double vc = d1 * d4 - d3 * d2;
    if(vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0)
    {
        const double v = d1 / (d1 - d3);
        const Eigen::Vector3d q = a + v * ab;
        return (p - q).squaredNorm();
    }

    const double vb = d5 * d2 - d1 * d6;
    if(vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0)
    {
        const double w = d2 / (d2 - d6);
        const Eigen::Vector3d q = a + w * ac;
        return (p - q).squaredNorm();
    }

    const double va = d3 * d6 - d5 * d4;
    if(va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
    {
        const double w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        const Eigen::Vector3d q = b + w * (c - b);
        return (p - q).squaredNorm();
    }

    const double denom = va + vb + vc;
    if(std::fabs(denom) < 1.0e-30)
    {
        return point_segment_sqdist(p, a, b);
    }
    const double inv = 1.0 / denom;
    const double v = vb * inv;
    const double w = vc * inv;
    const Eigen::Vector3d q = a + v * ab + w * ac;
    return (p - q).squaredNorm();
}

// Simple AABB for BVH nodes
struct BVH_AABB {
    Eigen::Vector3d min;
    Eigen::Vector3d max;
    
    BVH_AABB() : min(std::numeric_limits<double>::max(), 
                     std::numeric_limits<double>::max(), 
                     std::numeric_limits<double>::max()),
                 max(-std::numeric_limits<double>::max(), 
                     -std::numeric_limits<double>::max(), 
                     -std::numeric_limits<double>::max()) {}
    
    // Expand AABB to include a point
    void expand(const Eigen::Vector3d& point) {
        min.x() = std::min(min.x(), point.x());
        min.y() = std::min(min.y(), point.y());
        min.z() = std::min(min.z(), point.z());
        max.x() = std::max(max.x(), point.x());
        max.y() = std::max(max.y(), point.y());
        max.z() = std::max(max.z(), point.z());
    }
    
    // Expand to include another AABB
    void expand(const BVH_AABB& other) {
        min.x() = std::min(min.x(), other.min.x());
        min.y() = std::min(min.y(), other.min.y());
        min.z() = std::min(min.z(), other.min.z());
        max.x() = std::max(max.x(), other.max.x());
        max.y() = std::max(max.y(), other.max.y());
        max.z() = std::max(max.z(), other.max.z());
    }
    
    // Check if sphere intersects AABB
    bool intersects_sphere(const Eigen::Vector3d& center, double radius) const {
        double sqDist = 0.0;
        
        for(int i = 0; i < 3; ++i) {
            double v = center(i);
            if(v < min(i)) sqDist += (min(i) - v) * (min(i) - v);
            if(v > max(i)) sqDist += (v - max(i)) * (v - max(i));
        }
        
        return sqDist <= radius * radius;
    }
    
    // Get center of AABB
    Eigen::Vector3d center() const {
        return 0.5 * (min + max);
    }
    
    // Get surface area (SAH heuristic)
    double surface_area() const {
        Eigen::Vector3d d = max - min;
        return 2.0 * (d.x() * d.y() + d.y() * d.z() + d.z() * d.x());
    }
};

// Triangle structure for BVH
struct BVH_Triangle {
    Eigen::Vector3d v0, v1, v2;  // Vertices
    int triangle_id;              // Original triangle index
    
    BVH_Triangle() : triangle_id(-1) {}
    
    BVH_Triangle(const Eigen::Vector3d& a, const Eigen::Vector3d& b, 
                 const Eigen::Vector3d& c, int id) 
        : v0(a), v1(b), v2(c), triangle_id(id) {}
    
    // Get centroid of triangle
    Eigen::Vector3d centroid() const {
        return (v0 + v1 + v2) / 3.0;
    }
    
    // Get AABB of this triangle
    BVH_AABB get_aabb() const {
        BVH_AABB box;
        box.expand(v0);
        box.expand(v1);
        box.expand(v2);
        return box;
    }
};

// BVH Node
struct BVH_Node {
    BVH_AABB bbox;                      // Bounding box
    BVH_Node* left;                     // Left child
    BVH_Node* right;                    // Right child
    std::vector<int> triangle_indices;  // Triangle indices (only for leaf nodes)
    bool is_leaf;                       // Is this a leaf node?
    
    BVH_Node() : left(nullptr), right(nullptr), is_leaf(false) {}
    
    ~BVH_Node() {
        if(left) delete left;
        if(right) delete right;
    }
};

// BVH Tree class
class BVH_Tree {
private:
    BVH_Node* root;
    std::vector<BVH_Triangle> triangles;
    int max_triangles_per_leaf;
    
    // Recursive build function
    BVH_Node* build_recursive(std::vector<int>& triangle_indices, int depth) {
        BVH_Node* node = new BVH_Node();
        
        // Compute bounding box for this node
        for(int idx : triangle_indices) {
            node->bbox.expand(triangles[idx].get_aabb());
        }
        
        // Leaf node condition: few triangles or max depth
        if(triangle_indices.size() <= max_triangles_per_leaf || depth > 20) {
            node->is_leaf = true;
            node->triangle_indices = triangle_indices;
            return node;
        }
        
        // Find longest axis
        Eigen::Vector3d extent = node->bbox.max - node->bbox.min;
        int axis = 0;
        if(extent.y() > extent.x()) axis = 1;
        if(extent.z() > extent(axis)) axis = 2;
        
        // Sort triangles along this axis by their centroid
        std::sort(triangle_indices.begin(), triangle_indices.end(),
                  [this, axis](int a, int b) {
                      return triangles[a].centroid()(axis) < triangles[b].centroid()(axis);
                  });
        
        // Split at median
        int mid = triangle_indices.size() / 2;
        std::vector<int> left_indices(triangle_indices.begin(), triangle_indices.begin() + mid);
        std::vector<int> right_indices(triangle_indices.begin() + mid, triangle_indices.end());
        
        // Recursively build children
        node->left = build_recursive(left_indices, depth + 1);
        node->right = build_recursive(right_indices, depth + 1);
        
        return node;
    }
    
    // Recursive query function
    bool query_recursive(BVH_Node* node, const Eigen::Vector3d& center, double radius) const {
        if(!node) return false;
        
        // Test against node's bounding box
        if(!node->bbox.intersects_sphere(center, radius)) {
            return false;  // No intersection with this subtree
        }
        
        // Leaf: true sphere–triangle proximity (closest point on triangle)
        if(node->is_leaf) {
            const double r2 = radius * radius;
            for(int idx : node->triangle_indices) {
                const BVH_Triangle& tri = triangles[idx];
                const double d2 = point_triangle_sqdist(center, tri.v0, tri.v1, tri.v2);
                if(d2 <= r2)
                {
                    return true;
                }
            }
            return false;
        }
        
        // Internal node: recursively test children
        if(query_recursive(node->left, center, radius)) return true;
        if(query_recursive(node->right, center, radius)) return true;
        
        return false;
    }

public:
    BVH_Tree(int max_leaf_size = 4) 
        : root(nullptr), max_triangles_per_leaf(max_leaf_size) {}
    
    ~BVH_Tree() {
        if(root) delete root;
    }
    
    // Build BVH from triangle array
    void build(double** tri_x, double** tri_y, double** tri_z, int num_triangles) {
        if(root) {
            delete root;
            root = nullptr;
        }
        
        triangles.clear();
        triangles.reserve(num_triangles);
        
        // Load triangles
        for(int i = 0; i < num_triangles; ++i) {
            Eigen::Vector3d v0(tri_x[i][0], tri_y[i][0], tri_z[i][0]);
            Eigen::Vector3d v1(tri_x[i][1], tri_y[i][1], tri_z[i][1]);
            Eigen::Vector3d v2(tri_x[i][2], tri_y[i][2], tri_z[i][2]);
            triangles.emplace_back(v0, v1, v2, i);
        }
        
        // Build index list
        std::vector<int> all_indices(num_triangles);
        for(int i = 0; i < num_triangles; ++i) {
            all_indices[i] = i;
        }
        
        // Build tree
        root = build_recursive(all_indices, 0);
    }
    
    // Query: does sphere intersect any triangle?
    bool intersects_sphere(const Eigen::Vector3d& center, double radius) const {
        if(!root) return false;
        return query_recursive(root, center, radius);
    }
    
    // Get triangle count
    int triangle_count() const {
        return triangles.size();
    }
    
    // Check if built
    bool is_built() const {
        return root != nullptr;
    }
};

#endif // SIXDOF_BVH_H_


