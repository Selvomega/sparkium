#pragma once
#include "glm/glm.hpp"
#include "sparks/assets/vertex.h"
#include <vector>

namespace sparks {
struct AxisAlignedBoundingBox {
  float x_low{};
  float x_high{};
  float y_low{};
  float y_high{};
  float z_low{};
  float z_high{};
  std::vector<Vertex> vertices_; // The triangles (vertex triples) contained in the AABB. 
  AxisAlignedBoundingBox(float x_low,
                         float x_high,
                         float y_low,
                         float y_high,
                         float z_low,
                         float z_high);
  AxisAlignedBoundingBox(const glm::vec3 &position = glm::vec3{0.0f});
  AxisAlignedBoundingBox(const std::vector<Vertex> vertices);
  [[nodiscard]] bool IsIntersect(const glm::vec3 &origin,
                                 const glm::vec3 &direction,
                                 float t_min,
                                 float t_max) const;
  [[nodiscard]] float centroid(int dim) const;
  AxisAlignedBoundingBox operator&(const AxisAlignedBoundingBox &aabb) const;
  AxisAlignedBoundingBox operator|(const AxisAlignedBoundingBox &aabb) const;
  AxisAlignedBoundingBox &operator&=(const AxisAlignedBoundingBox &aabb);
  AxisAlignedBoundingBox &operator|=(const AxisAlignedBoundingBox &aabb);
};
}  // namespace sparks
