#include "sparks/assets/material.h"
#include <cmath>
#include <random>

#include "glm/fwd.hpp"
#include "glm/geometric.hpp"
#include "grassland/grassland.h"
#include "sparks/assets/scene.h"
#include "sparks/assets/texture.h"
#include "sparks/assets/util.h"
#include "sparks/util/util.h"

namespace sparks {

namespace {
std::unordered_map<std::string, MaterialType> material_name_map{
    {"lambertian", MATERIAL_TYPE_LAMBERTIAN},
    {"specular", MATERIAL_TYPE_SPECULAR},
    {"transmissive", MATERIAL_TYPE_TRANSMISSIVE},
    {"principled", MATERIAL_TYPE_PRINCIPLED},
    {"emission", MATERIAL_TYPE_EMISSION}};
}

Material::Material(Scene *scene, const tinyxml2::XMLElement *material_element)
    : Material() {
  if (!material_element) {
    return;
  }

  albedo_color = glm::vec3{1.0f};

  auto child_element = material_element->FirstChildElement("albedo");
  if (child_element) {
    albedo_color = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("albedo_texture");
  if (child_element) {
    std::string path = child_element->FindAttribute("value")->Value();
    Texture albedo_texture(1, 1);
    if (Texture::Load(path, albedo_texture)) {
      albedo_texture_id =
          scene->AddTexture(albedo_texture, PathToFilename(path));
    }
  }

  child_element = material_element->FirstChildElement("normal_texture");
  if (child_element) {
    std::string path = child_element->FindAttribute("value")->Value();
    Texture normal_texture(1, 1);
    if (Texture::Load(path, normal_texture)) {
      normal_texture_id =
          scene->AddTexture(normal_texture, PathToFilename(path));
    }
  }

  child_element = material_element->FirstChildElement("emission");
  if (child_element) {
    emission = StringToVec3(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("emission_strength");
  if (child_element) {
    emission_strength =
        std::stof(child_element->FindAttribute("value")->Value());
  }

  child_element = material_element->FirstChildElement("alpha");
  if (child_element) {
    alpha = std::stof(child_element->FindAttribute("value")->Value());
  }

  material_type =
      material_name_map[material_element->FindAttribute("type")->Value()];
}

Material::Material(const glm::vec3 &albedo) : Material() {
  albedo_color = albedo;
}

glm::vec3 Material::BRDF(const glm::vec3 &inDir, const glm::vec3 &outDir, const HitRecord &hit_record, const Scene* scene) const {
  /*
  Return the BRDF of the material. 
  scene is used to access the textures. 
  The `inDir` and `outDir` are the back-traced light directions. 
  */
  // TODO
  if (material_type==MATERIAL_TYPE_LAMBERTIAN) {
    // If the material is Lambertian type. 
    if (SameSideCheck(inDir, outDir, hit_record)) {
      return albedo_color * glm::vec3{scene->GetTextures()[albedo_texture_id].Sample(hit_record.tex_coord)} * float(1/M_PI);
    }
    return glm::vec3{0.0f};
  }
  else if (material_type==MATERIAL_TYPE_EMISSION) {
    // If the material is emissive.
    // Dealt like Lambertian material temporarily but without adding albedo color. 
    // `inDir` is not used here, so it is still fine even if invalid `hit_record.prev_direction` is passed in when explicitly sampling the light. 
    return emission * emission_strength;
  }
  else if (material_type==MATERIAL_TYPE_SPECULAR) {
    // If the material is specular.
    auto expected_out_dir = inDir - 2*glm::dot(inDir,hit_record.textured_normal)*hit_record.textured_normal;
    if (glm::dot(expected_out_dir,outDir)<1-(1e-3)) {
      // If the output direction is too different from the expected replected direction. 
      return glm::vec3{0.0f};
    }
    if (!SameSideCheck(inDir, outDir, hit_record)) {
      return glm::vec3{0.0f};
    }
    return glm::vec3{1.0f};
  }
  else if (material_type==MATERIAL_TYPE_TRANSMISSIVE) {
    return glm::vec3{0.0f};
  }
  else {
    // Other cases
    // TODO
    return glm::vec3{0.0f};
  }
}

}  // namespace sparks
