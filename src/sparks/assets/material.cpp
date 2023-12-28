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

#include <glm/glm.hpp>
#include <glm/gtc/random.hpp>

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
    if (SanityCheck(inDir, outDir, hit_record)) {
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
  else {
    // Other cases
    // TODO
    return glm::vec3{0.0f};
  }
}

std::pair<glm::vec3,float> Material::UniformSampling(const glm::vec3 &inDir, const HitRecord &hit_record) const {
  /*
  Sample a random direction in the hemisphere currently. 
  Do sanity check for 3 times. 
  Return the direction and the pdf value at the direction. 
  */
  // This is Uniform Sampling.
  glm::vec3 ret = glm::ballRand(1.0f);
  if (!SanityCheck(inDir, ret, hit_record)) {
    ret = -ret;
  }
  return std::make_pair(ret,1/(2*M_PI));
}

std::pair<glm::vec3, float> Material::ImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record) const {
  float exponent = 1.0f;
  // Cosine-weighted hemisphere sampling
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  float r1 = dis(gen); 
  float r2 = dis(gen);

  float phi = 2 * M_PI * r1;
  float x = cos(phi) * sqrt(r2);
  float y = sin(phi) * sqrt(r2);
  float z = sqrt(1 - r2);

  glm::vec3 reflectDir = glm::reflect(-inDir, hit_record.normal);

  // Transform the sample point to the coordinate system defined by the normal at hit_record
  glm::vec3 w = hit_record.normal;
  glm::vec3 u = glm::normalize(glm::cross((fabs(w.x) > 0.1f ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0)), w));
  glm::vec3 v = glm::cross(w, u);

  glm::vec3 sampleDir = glm::normalize(x * u + y * v + z * w);

  // Calculate the PDF for this direction
  // float pdf = glm::dot(reflectDir, sampleDir) / M_PI;
  float pdf = (exponent + 1) * pow(z, exponent) / (2 * M_PI);
  glm::vec3 ret = sampleDir;
  if (!SanityCheck(inDir, ret, hit_record)) {
    ret = -ret;
  }
  return std::make_pair(ret, pdf);
}




std::pair<glm::vec3, float> Material::CosImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record) const {
  // Cosine-weighted hemisphere sampling
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  float r1 = dis(gen); 
  float r2 = dis(gen);

  float phi = 2 * M_PI * r1;
  float x = cos(phi) * sqrt(r2);
  float y = sin(phi) * sqrt(r2);
  float z = sqrt(1 - r2);

  // Transform the sample point to the coordinate system defined by the normal at hit_record
  glm::vec3 w = hit_record.normal;
  glm::vec3 u = glm::normalize(glm::cross((fabs(w.x) > 0.1f ? glm::vec3(0, 1, 0) : glm::vec3(1, 0, 0)), w));
  glm::vec3 v = glm::cross(w, u);

  glm::vec3 sampleDir = glm::normalize(x * u + y * v + z * w);

  // Calculate the PDF for this direction
  float pdf = glm::dot(hit_record.normal, sampleDir) / M_PI;
  glm::vec3 ret = sampleDir;
  if (!SanityCheck(inDir, ret, hit_record)) {
    ret = -ret;
  }
  return std::make_pair(ret, pdf);
}
std::pair<glm::vec3, float> Material::MultiImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record) const {
  float cos_weight = 0.5f;
  float brdf_weight = 1 - cos_weight;
  auto cos_sample = Material::CosImportanceSampling(inDir, hit_record);
  auto brdf_sample = Material::ImportanceSampling(inDir, hit_record);
  // return std::make_pair(cos_sample.first*cos_weight+brdf_sample.first*brdf_weight, cos_sample.second*cos_weight+brdf_sample.second*brdf_weight);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 1.0);

  float r1 = dis(gen); 

  if (r1<cos_weight) {
    return cos_sample;
  }
  else {
    return brdf_sample;
  }
}



}  // namespace sparks