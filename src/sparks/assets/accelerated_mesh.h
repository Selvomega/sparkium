#pragma once
#include "sparks/assets/aabb.h"
#include "sparks/assets/mesh.h"

namespace sparks {

namespace {
struct TreeNode {
  AxisAlignedBoundingBox aabb{};
  int child[2]{-1, -1};
};
struct BVHNode {
  AxisAlignedBoundingBox aabb; // AABB
  std::unique_ptr<BVHNode> left; // left son
  std::unique_ptr<BVHNode> right; // right son
  int start; // the beginning index of the nodes under this node
  int end; // the end index of the nodes under this node
  BVHNode() : start(-1), end(-1) {}
};
}  // namespace

class AcceleratedMesh : public Mesh {
 public:
  AcceleratedMesh() = default;
  explicit AcceleratedMesh(const Mesh &mesh);
  AcceleratedMesh(const std::vector<Vertex> &vertices,
                  const std::vector<uint32_t> &indices);
  float TraceRay(const glm::vec3 &origin,
                 const glm::vec3 &direction,
                 float t_min,
                 HitRecord *hit_record) const override;
  void BuildAccelerationStructure();

 private:
  /*
   * You can add your acceleration structure contents here.
   * */
  // TODO 
};
}  // namespace sparks
