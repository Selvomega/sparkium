#pragma once
#include "cstdint"
#include "glm/glm.hpp"
#include "sparks/assets/util.h"
#include "sparks/assets/hit_record.h"

namespace sparks {

enum MaterialType : int {
  MATERIAL_TYPE_LAMBERTIAN = 0, // Uniformly reflect the lights
  MATERIAL_TYPE_SPECULAR = 1, // Total reflect
  MATERIAL_TYPE_TRANSMISSIVE = 2, // Lights can go through
  MATERIAL_TYPE_PRINCIPLED = 3,
  MATERIAL_TYPE_EMISSION = 4,
  MATERIAL_TYPE_METALLIC = 5,
  MATERIAL_TYPE_DIELECTRIC_GLOSSY = 6,
};

class Scene;

struct Material {
  glm::vec3 albedo_color{0.8f}; // The diffusive color of a material. 
  int albedo_texture_id{0};
  glm::vec3 emission{0.0f}; // The light emitted by the object. 
  float emission_strength{1.0f}; // The emission strength of the light emitted. 
  float alpha{1.0f}; // Measuring the opacity of the object. 
  MaterialType material_type{MATERIAL_TYPE_LAMBERTIAN}; // Material type. 
  float reserve[2]{}; // What is this used for?
  Material() = default;
  explicit Material(const glm::vec3 &albedo);
  Material(Scene *scene, const tinyxml2::XMLElement *material_element);

  [[nodiscard]] glm::vec3 BRDF(const glm::vec3 &inDir, const glm::vec3 &outDir, const HitRecord &hit_record, const Scene* scene) const;
  [[nodiscard]] glm::vec3 ctBRDF(const glm::vec3 &inDir, const glm::vec3 &outDir, const HitRecord &hit_record, const Scene* scene) const;
  [[nodiscard]] std::pair<glm::vec3,float> UniformSampling(const glm::vec3 &inDir, const HitRecord &hit_record) const;
  [[nodiscard]] std::pair<glm::vec3,float> ImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record) const;
  [[nodiscard]] std::pair<glm::vec3,float> CosImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record) const;
  [[nodiscard]] std::pair<glm::vec3,float> MultiImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record) const;
};
}  // namespace sparks
