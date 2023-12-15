#pragma once
#include "sparks/assets/aabb.h"
#include "sparks/assets/mesh.h"

#include <vector>

namespace sparks {

struct TreeNode {
  AxisAlignedBoundingBox aabb{};
  int child[2]{-1, -1};
};

class BVHTree {
  public:
    struct BVHNode {
      AxisAlignedBoundingBox aabb; // AABB
      std::vector<std::pair<const AxisAlignedBoundingBox *,int>> objs; // The list of objects in the BVH node. 
      std::shared_ptr<BVHNode> left; // left son
      std::shared_ptr<BVHNode> right; // right son
    };
    BVHTree();
    BVHTree(std::vector<AxisAlignedBoundingBox> &b);
    bool Intersect(const glm::vec3 &origin,
                   const glm::vec3 &direction,
                   float t_min,
                   std::vector<int> &result) const;
    AxisAlignedBoundingBox AccessBuf(int index) const;
  private:
    void Intersect_r(BVHNode *root,
                     const glm::vec3 &origin,
                     const glm::vec3 &direction,
                     float t_min,
                     std::vector<int> &result) const;
    BVHNode *Build(std::vector<std::pair<const AxisAlignedBoundingBox *,int>> &objs, int depth = 0);
    std::shared_ptr<BVHNode> root;
	  std::vector<AxisAlignedBoundingBox> buf;
};

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
  BVHTree tree; 
};
}  // namespace sparks
