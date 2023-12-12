#pragma once
#include "glm/glm.hpp"

namespace sparks {
struct HitRecord {
  int hit_entity_id{-1}; // The id of the entity hit
  glm::vec3 position{}; // The position of the intersection 
  glm::vec3 normal{}; // The normal after normal mapping or other processing
  glm::vec3 geometry_normal{}; // The normal in the model rather than "story"
  glm::vec3 tangent{}; 
  glm::vec2 tex_coord{}; // The 2D texture coordination 
  bool front_face{}; // Whether the hit point is at the front face of the object
};
}  // namespace sparks
