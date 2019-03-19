#ifndef RAY_TRACING_KD_TREE_H
#define RAY_TRACING_KD_TREE_H

#include <vector>
#include <memory>

#include "glm/geometric.hpp"
#include "triangles.h"

struct BoundingBox {
    glm::vec3 min;
    glm::vec3 max;
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
    Node(std::vector<int> triangles, int depth, BoundingBox BB, std::vector<Triangle> &triangles_) {
        depth_ = depth;
        int planeCoord = 0;
        BoundingBox leftBB = BB;
        BoundingBox rightBB = BB;
        float x = BB.max.x - BB.min.x;
        float y = BB.max.y - BB.min.y;
        float z = BB.max.z - BB.min.z;

        if (x > y && x > z) {
            plane_ = x / 2;
            planeCoord = 0;
        }

        if (y > x && y > z) {
            plane_ = y / 2;
            planeCoord = 1;
        }

        if (z > x && z > y) {
            plane_ = z / 2;
            planeCoord = 2;
        }

        leftBB.max[planeCoord] = plane_;
        rightBB.min[planeCoord] = plane_;

        std::vector<int> leftTriangles;
        std::vector<int> rightTriangles;

        for (int i: triangles) {
            int minCoord = std::min(triangles_[i].GetVertices()[0][planeCoord], triangles_[i].GetVertices()[1][planeCoord],
                             triangles_[i].GetVertices()[2][planeCoord]);
            int maxCoord = std::max(triangles_[i].GetVertices()[0][planeCoord], triangles_[i].GetVertices()[1][planeCoord],
                                    triangles_[i].GetVertices()[2][planeCoord]);
            if (minCoord <= plane_) {
                leftTriangles.push_back(i);
            }

            if (maxCoord >= plane_) {
                rightTriangles.push_back(i);
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
