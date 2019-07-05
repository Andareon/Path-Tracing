#ifndef RAY_TRACING_KD_TREE_H
#define RAY_TRACING_KD_TREE_H

#include <memory>

#include "glm/geometric.hpp"
#include "Tracer.h"

struct BoundingBox {
    glm::vec3 min;
    glm::vec3 max;

    BoundingBox() : max(std::numeric_limits<float>::min()), min(std::numeric_limits<float>::max()) {}

    explicit BoundingBox(const Triangle& triangle) :
        max(-std::numeric_limits<float>::max()),
        min(std::numeric_limits<float>::max())
    {
        for (const auto& vertex : triangle.GetVertices()) {
            for (int i = 0; i < 3; ++i) {
                max[i] = std::fmax(vertex[i], max[i]);
                min[i] = std::fmin(vertex[i], min[i]);
            }
        }
    }

    BoundingBox(glm::vec3 mn, glm::vec3 mx): min(mn), max(mx){}
};

class BaseNode {
protected:
    BoundingBox BB_;
public:
    virtual bool Trace(Ray &ray, IntersectionOptions &options, float dist) = 0;
    virtual bool isEmpty() = 0;
    virtual ~BaseNode() = default;
};

class Leave : public BaseNode {
private:
    std::vector<int> trianglesIndex_;
    std::vector<Triangle> &triangles_;
public:
    Leave(std::vector<int> trianglesIndex, std::vector<Triangle> &triangles) : trianglesIndex_(std::move(trianglesIndex)),
                                                                               triangles_(triangles){}
    bool Trace(Ray &ray, IntersectionOptions &options, float dist) {
        float distance = dist;
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

    bool isEmpty() {
        return trianglesIndex_.empty();
    }
};

class Node : public BaseNode {
private:
    std::unique_ptr<BaseNode> left_;
    std::unique_ptr<BaseNode> right_;
    int planeCoord_;
    float plane_;
public:
    Node(const std::vector<int> &trianglesIndex_, int depth, BoundingBox BB, std::vector<Triangle> &triangles_) {
        BB_ = BB;
        BoundingBox leftBB = BB;
        BoundingBox rightBB = BB;
        const glm::vec3 bbSides = BB.max - BB.min;
        if (bbSides.x >= bbSides.y && bbSides.x >= bbSides.z) {
            planeCoord_ = 0;
        } else if (bbSides.y >= bbSides.x && bbSides.y >= bbSides.z) {
            planeCoord_ = 1;
        } else {
            planeCoord_ = 2;
        }

        const int SPACE_PARTITIONS = 21;
        float currentPlane = BB.min[planeCoord_];
        float delta = bbSides[planeCoord_] / SPACE_PARTITIONS;
        plane_ = 0;
        float minSAH = trianglesIndex_.size() * (SPACE_PARTITIONS + 1);
        bool splited = false;
        for (int leftPartitions = 1; leftPartitions < SPACE_PARTITIONS; ++leftPartitions) {
            currentPlane += delta;
            int leftCount = 0;
            int rightCount = 0;
            for (int j = 0; j < trianglesIndex_.size(); ++j) {
                if (BoundingBox(triangles_[j]).min[planeCoord_] < currentPlane) {
                    leftCount++;
                }
                if (BoundingBox(triangles_[j]).max[planeCoord_] > currentPlane) {
                    rightCount++;
                }
            }
            float SAH = leftCount * leftPartitions + rightCount * (SPACE_PARTITIONS - leftPartitions);
            if (SAH < minSAH) {
                minSAH = SAH;
                plane_ = currentPlane;
                splited = true;
            }
        }
        if (!splited) {
            plane_ = BB.min[planeCoord_] + bbSides[planeCoord_] / 2;
        }


        leftBB.max[planeCoord_] = plane_;
        rightBB.min[planeCoord_] = plane_;

        std::vector<int> leftTriangles;
        std::vector<int> rightTriangles;

        for (int triangle_id: trianglesIndex_) {
            const BoundingBox triangleBB = BoundingBox(triangles_[triangle_id]);
            if (triangleBB.min[planeCoord_] <= plane_) {
                leftTriangles.push_back(triangle_id);
            }

            if (triangleBB.max[planeCoord_] >= plane_) {
                rightTriangles.push_back(triangle_id);
            }
        }

        if (!splited || leftTriangles.size() == trianglesIndex_.size() || rightTriangles.size() == trianglesIndex_.size()) {
            left_ = std::make_unique<Leave>(trianglesIndex_, triangles_);
            right_ = std::make_unique<Leave>(std::vector<int>(), triangles_);
            return;
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

    bool isEmpty() {
        return false;
    }

    bool Trace(Ray &ray, IntersectionOptions &options, float dist) {
        float t = (plane_ - ray.GetBegin()[planeCoord_]) / ray.GetDirection()[planeCoord_];
        float tNear = -INFINITY, tFar = INFINITY;
        for (int i = 0; i < 3; ++i) {
            if (fabs(ray.GetDirection()[i]) > 1e-7) {
                float t1 = (BB_.min[i] - ray.GetBegin()[i]) / ray.GetDirection()[i];
                float t2 = (BB_.max[i] - ray.GetBegin()[i]) / ray.GetDirection()[i];
                if (t1 > t2) {
                    std::swap(t1, t2);
                }
                if (t1 > tNear) {
                    tNear = t1;
                }
                if (t2 < tFar) {
                    tFar = t2;
                }
            }
        }
        if (tNear < 0) {
            tNear = 0;
        }
        if (tNear >= tFar) {
            return false;
        }

        if (right_->isEmpty()) {
            return left_->Trace(ray, options, tFar);
        }
        if (t > tFar) {
            if (ray.GetDirection()[planeCoord_] > 0) {
                return left_->Trace(ray, options, tFar);
            }
            return right_->Trace(ray, options, tFar);
        } else if (t < tNear) {
            if (ray.GetDirection()[planeCoord_] > 0) {
                return right_->Trace(ray, options, tFar);
            }
            return left_->Trace(ray, options, tFar);
        } else {
            if (ray.GetBegin()[planeCoord_] < plane_) {
                return left_->Trace(ray, options, tFar) || right_->Trace(ray, options, tFar);
            } else if (ray.GetBegin()[planeCoord_] > plane_){
                return right_->Trace(ray, options, tFar) || left_->Trace(ray, options, tFar);
            } else {
                if (ray.GetDirection()[planeCoord_] > 0) {
                    return right_->Trace(ray, options, tFar);// || left_->Trace(ray, options);
                } else {
                    return left_->Trace(ray, options, tFar);// || right_->Trace(ray, options);
                }
            }
        }
    }
};

class KDTreeTracer : public BaseTracer {
    Node root;
public:
    KDTreeTracer(const std::vector<int> &trianglesIndex_, int depth, BoundingBox BB, std::vector<Triangle> &triangles_) :
    root(trianglesIndex_, depth, BB, triangles_) {}

    bool Trace(Ray &ray, IntersectionOptions &options) final {
        return root.Trace(ray, options, INFINITY);
    }
};

#endif //RAY_TRACING_KD_TREE_H
