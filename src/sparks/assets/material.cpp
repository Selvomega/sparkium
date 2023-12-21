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
      return albedo_color * glm::vec3{scene->GetTextures()[albedo_texture_id].Sample(hit_record.tex_coord)} * float(1/3.14);
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
  // TODO
  glm::vec3 ret = glm::ballRand(1.0f);
  if (!SanityCheck(inDir, ret, hit_record)) {
   ret = -ret;
  }
return std::make_pair(ret,1/(2*3.14));
}


//std::pair<glm::vec3, float> Material::ImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record) const {    float phongExponent = 100; // ���ݲ�����������    glm::vec3 reflectDir = glm::reflect(-inDir, hit_record.normal);    reflectDir = glm::normalize(reflectDir);    glm::vec3 w_vec = reflectDir;    glm::vec3 u_vec = hit_record.tangent;    glm::vec3 v_vec = glm::cross(w_vec, u_vec);    glm::vec3 sampledDir;    float pdf;    for (int i = 0; i < 3; i++) {        // ʹ��Phong�ֲ�ģ�Ͳ���        std::random_device rd;        std::mt19937 gen(rd());        std::uniform_real_distribution<> dis(0, 1);        double u = dis(gen);        double v = dis(gen);        double theta = std::acos(std::pow(u, 1 / (phongExponent + 1)));        double phi = 2 * M_PI * v;        // ת�����ֲ�����ϵ        float x = std::sin(theta) * std::cos(phi);        float y = std::sin(theta) * std::sin(phi);        float z = std::cos(theta);        sampledDir = glm::normalize(x * u_vec + y * v_vec + z * w_vec);        if (SanityCheck(sampledDir, hit_record)) {            float cosine = std::max(glm::dot(sampledDir, reflectDir), 0.0f);            pdf = (phongExponent + 1) / (2 * M_PI) * std::pow(cosine, phongExponent);            return std::make_pair(sampledDir, pdf);        }    }    // try 3 times at most    float cosine = std::max(glm::dot(sampledDir, reflectDir), 0.0f);    pdf = (phongExponent + 1) / (2 * M_PI) * std::pow(cosine, phongExponent);    return std::make_pair(sampledDir, pdf);}


//std::pair<glm::vec3, float> Material::UniformSampling(
//    const glm::vec3 &inDir,
//    const HitRecord &hit_record) const {
//    glm::vec3 ret{};
//  for (int i = 0; i< 3; i++) {
//    // Sample uniformly over the hemisphere.
//    std::random_device rd;
//    std::mt19937 gen(rd());
//    std::uniform_real_distribution<> dis(0, 1);
//    double u = dis(gen);
//    double v = dis(gen);
//    double theta = acos(sqrt(1.0 - u));
//    double phi = 2 *3.14 * v;
//    float x = sin(theta) * cos(phi);
//    float y = sin(theta) * sin(phi);
//    float z = cos(theta);
//    glm::vec3 w_vec = hit_record.normal;
//    glm::vec3 u_vec = hit_record.tangent;
//    glm::vec3 v_vec = glm::cross(w_vec, u_vec);
//    ret = glm::normalize(x* u_vec + y* v_vec + z* w_vec);
//    if (SanityCheck(inDir, ret, hit_record)) {
//      return std::make_pair(ret, 1 / (2* 3.14));
//    }
//  }
//  return std::make_pair(ret, 1 / (2* 3.14));
//}


std::pair<glm::vec3, float> Material::ImportanceSampling(
    const glm::vec3 &inDir, const HitRecord &hit_record) const {
    float phongExponent = 100;  // ���ݲ�����������
    glm::vec3 reflectDir = glm::reflect(-inDir, hit_record.normal);
    reflectDir = glm::normalize(reflectDir);
    glm::vec3 w_vec = reflectDir;
    glm::vec3 u_vec = hit_record.tangent;
    glm::vec3 v_vec = glm::cross(w_vec, u_vec);
    glm::vec3 sampledDir;
    float pdf;
    for (int i = 0; i < 3;
         i++) {  // ʹ��Phong�ֲ�ģ�Ͳ���       
        std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_real_distribution<> dis(0, 1);
      double u = dis(gen);
      double v = dis(gen);
      double theta = std::acos(std::pow(u, 1 / (phongExponent + 1)));
      double phi = 2 * 3.14 * v;  // ת�����ֲ�����ϵ       
      float x = std::sin(theta) * std::cos(phi);
      float y = std::sin(theta) * std::sin(phi);
      float z = std::cos(theta);
      sampledDir = glm::normalize(x * u_vec + y * v_vec + z * w_vec);
      if (SanityCheck(inDir, sampledDir, hit_record)) {
        float cosine = std::max(glm::dot(sampledDir, reflectDir), 0.0f);
        pdf =
            (phongExponent + 1) / (2 * 3.14) * std::pow(cosine, phongExponent);
        return std::make_pair(sampledDir, pdf);
      }
    }  // try 3 times at most
    float cosine = std::max(glm::dot(sampledDir, reflectDir), 0.0f);
    pdf = (phongExponent + 1) / (2 * 3.14) * std::pow(cosine, phongExponent);
    return std::make_pair(sampledDir, pdf);
  }

std::pair<glm::vec3, float> Material::CosImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record) const {
  float theta, phi, rand1, rand2;
  rand1 = fract(RandomFloat());
  rand2 = fract(RandomFloat());

  theta = acos(sqrt(1-rand1));
  phi = 2*PI*rand2;

  vec3 out_direction;

  vec3 local_z = normalize(hit_record.normal);
  vec3 local_x = normalize(direction - dot(direction, local_z) * local_z);
  vec3 local_y = cross(local_z,local_x);

  out_direction = normalize(sin(theta)*cos(phi)*local_x + sin(theta)*sin(phi)*local_y + cos(theta)*local_z);
  return std::make_pair(out_direction, cos(theta)/3.14);
}
std::pair<glm::vec3, float> Material::MultiImportanceSampling(
    const glm::vec3 &inDir,
    const HitRecord &hit_record) const { 
    // TODO
  float phongExponent = 100;
  float weightPhong = 0.5;
  float weightUniform = 0.5;
  glm::vec3 sampleDirPhong, sampleDirUniform;
  float pdfPhong, pdfUniform;
  std::tie(sampleDirPhong, pdfPhong) = ImportanceSampling(inDir, hit_record);
  std::tie(sampleDirUniform, pdfUniform) = UniformSampling(inDir, hit_record);
  float weightSum = weightPhong * pdfPhong + weightUniform * pdfUniform;
  glm::vec3 combinedDir = (weightPhong * pdfPhong * sampleDirPhong +
                           weightUniform * pdfUniform * sampleDirUniform) /
                          weightSum;
  float combinedPdf =
      (pdfPhong * weightPhong + pdfUniform * weightUniform) / weightSum;
  return std::make_pair(combinedDir, combinedPdf);
}

}  // namespace sparks
