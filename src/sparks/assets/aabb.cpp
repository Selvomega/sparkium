#include "sparks/assets/aabb.h"
#include <cmath>

#include "algorithm"
#include "glm/fwd.hpp"
#include "grassland/grassland.h"
#include "random"

namespace sparks {

AxisAlignedBoundingBox::AxisAlignedBoundingBox(float x_low,
                                               float x_high,
                                               float y_low,
                                               float y_high,
                                               float z_low,
                                               float z_high) {
  this->x_low = x_low;
  this->x_high = x_high;
  this->y_low = y_low;
  this->y_high = y_high;
  this->z_low = z_low;
  this->z_high = z_high;
}

AxisAlignedBoundingBox::AxisAlignedBoundingBox(const glm::vec3 &position) {
  x_low = position.x;
  x_high = position.x;
  y_low = position.y;
  y_high = position.y;
  z_low = position.z;
  z_high = position.z;
}

AxisAlignedBoundingBox::AxisAlignedBoundingBox(const std::vector<Vertex> vertices) {
  assert(vertices.size()%3==0);
  assert(vertices.size()!=0);
  x_low = vertices[0].position.x;
  x_high = vertices[0].position.x;
  y_low = vertices[0].position.y;
  y_high = vertices[0].position.y;
  z_low = vertices[0].position.z;
  z_high = vertices[0].position.z;
  vertices_.reserve(vertices.size());
  for (int i=0; i<vertices.size(); i++) {
    vertices_.emplace_back(vertices[i]); // Put the vertices into the 
    x_low = std::min(x_low,vertices[i].position.x);
    x_high = std::max(x_high,vertices[i].position.x);
    y_low = std::min(y_low,vertices[i].position.y);
    y_high = std::max(y_high,vertices[i].position.y);
    z_low = std::min(z_low,vertices[i].position.z);
    z_high = std::max(z_high,vertices[i].position.z);
  }
}

bool AxisAlignedBoundingBox::IsIntersect(const glm::vec3 &origin,
                                         const glm::vec3 &direction,
                                         float t_min,
                                         float t_max) const {
  /*
  Check whether the ray with given range of parameters has a intersection with the bounding box.
  */
  if (x_low <= origin.x && origin.x <= x_high && y_low <= origin.y &&
      origin.y <= y_high && z_low <= origin.z && origin.z <= z_high) {
    return true;
  }
  float intersection_range_low = t_max * (1.0f + t_min);
  float intersection_range_high = 0.0f;
  float t;
  glm::vec3 intersection;
#define TestIntersection(x, y, z)                                     \
  if (std::abs(direction.x) > 1e-5) {                                 \
    float inv_d = 1.0f / direction.x;                                 \
    t = (x##_low - origin.x) * inv_d;                                 \
    intersection = origin + direction * t;                            \
    if (y##_low <= intersection.y && intersection.y <= y##_high &&    \
        z##_low <= intersection.z && intersection.z <= z##_high) {    \
      intersection_range_low = std::min(intersection_range_low, t);   \
      intersection_range_high = std::max(intersection_range_high, t); \
    }                                                                 \
    t = (x##_high - origin.x) * inv_d;                                \
    intersection = origin + direction * t;                            \
    if (y##_low <= intersection.y && intersection.y <= y##_high &&    \
        z##_low <= intersection.z && intersection.z <= z##_high) {    \
      intersection_range_low = std::min(intersection_range_low, t);   \
      intersection_range_high = std::max(intersection_range_high, t); \
    }                                                                 \
  }
  TestIntersection(x, y, z);
  TestIntersection(z, x, y);
  TestIntersection(y, z, x);
  return intersection_range_high >= t_min && intersection_range_low <= t_max;
}

float AxisAlignedBoundingBox::centroid(int dim) const {
  /*
  Return the centeroid of the AABB. 
  */
  if (dim==0) {
    return (x_low+x_high)/2;
  }
  else if (dim==1) {
    return (y_low+y_high)/2;
  }
  return (z_low+z_high)/2;
}

AxisAlignedBoundingBox AxisAlignedBoundingBox::operator&(
    const AxisAlignedBoundingBox &aabb) const {
  return {std::max(x_low, aabb.x_low), std::min(x_high, aabb.x_high),
          std::max(y_low, aabb.y_low), std::min(y_high, aabb.y_high),
          std::max(z_low, aabb.z_low), std::min(z_high, aabb.z_high)};
}

AxisAlignedBoundingBox AxisAlignedBoundingBox::operator|(
    const AxisAlignedBoundingBox &aabb) const {
  return {std::min(x_low, aabb.x_low), std::max(x_high, aabb.x_high),
          std::min(y_low, aabb.y_low), std::max(y_high, aabb.y_high),
          std::min(z_low, aabb.z_low), std::max(z_high, aabb.z_high)};
}

AxisAlignedBoundingBox &AxisAlignedBoundingBox::operator&=(
    const AxisAlignedBoundingBox &aabb) {
  x_low = std::max(x_low, aabb.x_low);
  x_high = std::min(x_high, aabb.x_high);
  y_low = std::max(y_low, aabb.y_low);
  y_high = std::min(y_high, aabb.y_high);
  z_low = std::max(z_low, aabb.z_low);
  z_high = std::min(z_high, aabb.z_high);
  return *this;
}

AxisAlignedBoundingBox &AxisAlignedBoundingBox::operator|=(
    const AxisAlignedBoundingBox &aabb) {
  x_low = std::min(x_low, aabb.x_low);
  x_high = std::max(x_high, aabb.x_high);
  y_low = std::min(y_low, aabb.y_low);
  y_high = std::max(y_high, aabb.y_high);
  z_low = std::min(z_low, aabb.z_low);
  z_high = std::max(z_high, aabb.z_high);
  return *this;
}

glm::vec3 AxisAlignedBoundingBox::Rand() const {
  /*
  Randomly choose a point in the bounding box. 
  */
  std::random_device nrd;
  std::mt19937 gen(nrd());
  std::uniform_real_distribution<> dis(0.0, 1.0);
  float x_rand = dis(gen);
  float y_rand = dis(gen);
  float z_rand = dis(gen);
  return glm::vec3{(x_high-x_low)*x_rand+x_low,
                   (y_high-y_low)*y_rand+y_low,
                   (z_high-z_low)*z_rand+z_low};
}

}  // namespace sparks
