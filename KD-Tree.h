#ifndef RAY_TRACING_KD_TREE_H
#define RAY_TRACING_KD_TREE_H

#include <vector>
#include <memory>

#include "glm/geometric.hpp"
#include "triangles.h"

struct BoundingBox {
    glm::vec3 min;
    glm::vec3 max;

    BoundingBox() : max(std::numeric_limits<float>::min()), min(std::numeric_limits<float>::max()) {}
    explicit BoundingBox(const Triangle& triangle) :
        max(std::numeric_limits<float>::min()),
        min(std::numeric_limits<float>::max())
    {
        for (const auto& vertex : triangle.GetVertices()) {
            max = glm::max(vertex, max);
            min = glm::min(vertex, min);
        }
    }
};

class BaseNode {
protected:
    BoundingBox BB_;
public:
    virtual ~BaseNode() = default;
};

class Leave : public BaseNode {
private:
    std::vector<int> trianglesIndex_;
    std::vector<Triangle> &triangles_;
public:
    Leave(std::vector<int> trianglesIndex, std::vector<Triangle> &triangles) : trianglesIndex_(std::move(trianglesIndex)),
    triangles_(triangles){}
};

class Node : public BaseNode {
private:
    std::unique_ptr<BaseNode> left_;
    std::unique_ptr<BaseNode> right_;
    int depth_;
    int planeCoord_;
    float plane_;
public:
    Node(const std::vector<int> &triangles, int depth, BoundingBox BB, std::vector<Triangle> &triangles_) {
        depth_ = depth;
        int planeCoord = 0;
        BoundingBox leftBB = BB;
        BoundingBox rightBB = BB;
        const glm::vec3 bbSides = BB.max - BB.min;
        if (bbSides.x > bbSides.y && bbSides.x > bbSides.z) {
            planeCoord = 0;
        } else if (bbSides.y > bbSides.x && bbSides.y > bbSides.z) {
            planeCoord = 1;
        } else {
            planeCoord = 2;
        }

        plane_ = bbSides[planeCoord] * 0.5;
        leftBB.max[planeCoord] = plane_;
        rightBB.min[planeCoord] = plane_;

        std::vector<int> leftTriangles;
        std::vector<int> rightTriangles;

        for (int triangle_id: triangles) {
            const BoundingBox triangleBB = BoundingBox(triangles_[triangle_id]);
            if (triangleBB.min[planeCoord] <= plane_) {
                leftTriangles.push_back(triangle_id);
            }

            if (triangleBB.max[planeCoord] >= plane_) {
                rightTriangles.push_back(triangle_id);
            }
        }

        if (depth == 1 || leftTriangles.size() <= 1) {
            left_ = std::make_unique<Leave>(leftTriangles, triangles);
        } else {
            left_ = std::make_unique<Node>(leftTriangles, depth - 1, leftBB, triangles_);
        }

        if (depth == 1 || rightTriangles.size() <= 1) {
            right_ = std::make_unique<Leave>(rightTriangles, triangles);
        } else {
            right_ = std::make_unique<Node>(rightTriangles, depth - 1, rightBB, triangles_);
        }
    }
};
#endif //RAY_TRACING_KD_TREE_H
