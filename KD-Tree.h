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

    BoundingBox(glm::vec3 mn, glm::vec3 mx): min(mn), max(mx){}
};

struct IntersectionOptions {
    glm::vec4 intersectPoint_;
    glm::vec4 normal_;
    int materialIndex_;
};

class BaseNode {
protected:
    BoundingBox BB_;
public:
    virtual bool Trace(Ray &ray, IntersectionOptions &options) = 0;
    virtual ~BaseNode() = default;
};

class Leave : public BaseNode {
private:
    std::vector<int> trianglesIndex_;
    std::vector<Triangle> &triangles_;
public:
    Leave(std::vector<int> trianglesIndex, std::vector<Triangle> &triangles) : trianglesIndex_(std::move(trianglesIndex)),
                                                                               triangles_(triangles){}
    bool Trace(Ray &ray, IntersectionOptions &options) {
        float distance = INFINITY;
        int nearest = -1;
        for (int i: trianglesIndex_) {
            if (triangles_[i].Intersect(ray, distance)) {
                nearest = i;
            }
        }
        if (nearest != -1) {
            options.intersectPoint_ = ray.GetBegin() + ray.GetDirection() * distance;
            options.materialIndex_ = triangles_[nearest].GetMaterial();
            options.normal_ = triangles_[nearest].GetNormal();
            return true;
        }
        return false;
    }
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
        BoundingBox leftBB = BB;
        BoundingBox rightBB = BB;
        const glm::vec3 bbSides = BB.max - BB.min;
        if (bbSides.x > bbSides.y && bbSides.x > bbSides.z) {
            planeCoord_ = 0;
        } else if (bbSides.y > bbSides.x && bbSides.y > bbSides.z) {
            planeCoord_ = 1;
        } else {
            planeCoord_ = 2;
        }

        plane_ = bbSides[planeCoord_] * 0.5;
        leftBB.max[planeCoord_] = plane_;
        rightBB.min[planeCoord_] = plane_;

        std::vector<int> leftTriangles;
        std::vector<int> rightTriangles;

        for (int triangle_id: triangles) {
            const BoundingBox triangleBB = BoundingBox(triangles_[triangle_id]);
            if (triangleBB.min[planeCoord_] <= plane_) {
                leftTriangles.push_back(triangle_id);
            }

            if (triangleBB.max[planeCoord_] >= plane_) {
                rightTriangles.push_back(triangle_id);
            }
        }

        if (depth == 1 || leftTriangles.size() <= 1) {
            left_ = std::make_unique<Leave>(leftTriangles, triangles_);
        } else {
            left_ = std::make_unique<Node>(leftTriangles, depth - 1, leftBB, triangles_);
        }

        if (depth == 1 || rightTriangles.size() <= 1) {
            right_ = std::make_unique<Leave>(rightTriangles, triangles_);
        } else {
            right_ = std::make_unique<Node>(rightTriangles, depth - 1, rightBB, triangles_);
        }
    }

    bool Trace(Ray &ray, IntersectionOptions &options) {
        float t = (plane_ - ray.GetBegin()[planeCoord_]) / ray.GetDirection()[planeCoord_];
        float tNear = (BB_.min[planeCoord_] - ray.GetBegin()[planeCoord_]) / ray.GetDirection()[planeCoord_];
        float tFar = (BB_.max[planeCoord_] - ray.GetBegin()[planeCoord_]) / ray.GetDirection()[planeCoord_];
        if (tNear > tFar) {
            std::swap(tNear, tFar);
        }
        if (t >= tFar) {
            return left_->Trace(ray, options);
        } else if (t <= tNear) {
            return right_->Trace(ray, options);
        } else {
            if (ray.GetBegin()[planeCoord_] < plane_) {
                return left_->Trace(ray, options) || right_->Trace(ray, options);
            } else if (ray.GetBegin()[planeCoord_] > plane_){
                return right_->Trace(ray, options) || left_->Trace(ray, options);
            } else {
                if (ray.GetDirection()[planeCoord_] > 0) {
                    return right_->Trace(ray, options) || left_->Trace(ray, options);

                } else {
                    return left_->Trace(ray, options) || right_->Trace(ray, options);

                }
            }
        }
    }
};
#endif //RAY_TRACING_KD_TREE_H
