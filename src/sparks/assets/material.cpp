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
    {"emission", MATERIAL_TYPE_EMISSION},
    {"metallic", MATERIAL_TYPE_METALLIC},
    {"glossy", MATERIAL_TYPE_DIELECTRIC_GLOSSY},
    };
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



glm::vec3 Material::ctBRDF(const glm::vec3 &inDir, const glm::vec3 &outDir, const HitRecord &hit_record, const Scene* scene) const {
  /*
  Return the BRDF of the material. 
  scene is used to access the textures. 
  The `inDir` and `outDir` are the back-traced light directions. 
  */
  // Cook-Torrance BRDF
  // TODO
  if (material_type==MATERIAL_TYPE_LAMBERTIAN) {
    // If the material is Lambertian type. 
    // TODO: Lambertian BRDF
    if (SameSideCheck(inDir, outDir, hit_record)) {
      return albedo_color * glm::vec3{scene->GetTextures()[albedo_texture_id].Sample(hit_record.tex_coord)} * float(1/M_PI);
    }
    return glm::vec3{0.0f};
  }
  else if (material_type==MATERIAL_TYPE_EMISSION) {
    // If the material is emissive.
    // Dealt like Lambertian material temporarily but without adding albedo color. 
    // `inDir` is not used here, so it is still fine even if invalid `hit_record.prev_direction` is passed in when explicitly sampling the light. 
    // TODO : Emission BRDF
    return emission * emission_strength;
  }
  else if (material_type == MATERIAL_TYPE_METALLIC) {
    // Cook-Torrance BRDF for metallic surfaces
    if (SameSideCheck(inDir, outDir, hit_record)) {
      glm::vec3 normal = hit_record.textured_normal;
      glm::vec3 halfway = glm::normalize(-inDir + outDir);

      // Specular reflection (Fresnel-Schlick approximation)
      float cosTheta = glm::clamp(glm::dot(normal, halfway), 0.0f, 1.0f);
      // float cosTheta = glm::clamp(-glm::dot(normal, inDir), 0.0f, 1.0f);
      // glm::vec3 F0 = glm::vec3(0.04f); // Assuming non-metallic surface
      glm::vec3 F0 = albedo_color; 
      glm::vec3 fresnel = F0 + (1.0f - F0) * pow(1.0f - cosTheta, 5.0f);
      // fresnel become (1,1,1), why?
      // Microfacet Distribution (Trowbridge-Reitz GGX) 
      float roughness = 0.1; // Roughness value of the surface
      float alpha = roughness * roughness;
      float NdotH = glm::max(glm::dot(normal, halfway), 0.0f);
      float D = alpha * alpha / (M_PI * pow(NdotH * NdotH * (alpha * alpha - 1.0f) + 1.0f, 2));
      // float D = 1; // GGX Distribution function

      // Geometric Attenuation (Schlick-GGX approximation)
      float k = (roughness + 1.0f) * (roughness + 1.0f) / 8.0f;
      float NdotL = glm::max(glm::dot(normal, inDir), 0.0f);
      float NdotV = glm::max(glm::dot(normal, outDir), 0.0f);
      NdotL = glm::max(NdotL, -glm::dot(normal, inDir));
      NdotV = glm::max(NdotV, -glm::dot(normal, outDir));

      float G1 = NdotL / (NdotL * (1.0f - k) + k);
      float G2 = NdotV / (NdotV * (1.0f - k) + k);
      float G = G1 * G2;
      // Final BRDF
      glm::vec3 numerator = fresnel * D * G;
      float denominator = 4 * glm::max(glm::dot(normal, inDir), 0.0f) * glm::max(glm::dot(normal, outDir), 0.0f) + 0.001f; // Adding a small value to avoid division by zero

      glm::vec3 specular = numerator / denominator;

      // Combine diffuse and specular
      float max_specular = glm::max(specular[0], glm::max(specular[1], specular[2]));
      // specular = specular / max_specular;
      specular = specular * 2.0f;
      return specular; 
    }
    return glm::vec3{0.0f};
  }
  else if (material_type == MATERIAL_TYPE_DIELECTRIC_GLOSSY) {
    // If the material is a glossy dielectric (non-metallic with specular highlights).
    if (SameSideCheck(inDir, outDir, hit_record)) {
        // Albedo color can be either a constant or obtained from a texture
        glm::vec3 albedo = albedo_color; // Assume albedo_color is a glm::vec3 representing the base color
        if (albedo_texture_id >= 0) {
            // If there is a texture, use it to modulate the albedo
            albedo *= glm::vec3{scene->GetTextures()[albedo_texture_id].Sample(hit_record.tex_coord)};
        }

        // Lambertian (Diffuse) component
        glm::vec3 diffuse = albedo * float(1.0 / M_PI);
        glm::vec3 normal = hit_record.textured_normal;
        // Specular (Cook-Torrance) component
        glm::vec3 halfway = glm::normalize(inDir + outDir);
        float roughness = 0.2; // Roughness value of the surface
        float alpha = roughness * roughness;
        float NdotH = glm::max(glm::dot(normal, halfway), 0.0f);
        float D = alpha * alpha / (M_PI * pow(NdotH * NdotH * (alpha * alpha - 1.0f) + 1.0f, 2));

        float k = (roughness + 1.0f) * (roughness + 1.0f) / 8.0f;
        float NdotL = glm::max(glm::dot(normal, inDir), 0.0f);
        float NdotV = glm::max(glm::dot(normal, outDir), 0.0f);
        NdotV = glm::max(NdotV, -glm::dot(normal, outDir));
        float G1 = NdotL / (NdotL * (1.0f - k) + k);
        float G2 = NdotV / (NdotV * (1.0f - k) + k);
        float G = G1 * G2;
        glm::vec3 F0 = glm::vec3(0.04f); // Fresnel at normal incidence for non-metals; can be overridden based on material properties
        float cosTheta = glm::clamp(glm::dot(halfway, outDir), 0.0f, 1.0f);
        glm::vec3 F = F0 + (glm::vec3(1.0f) - F0) * pow(1.0f - cosTheta, 5.0f);

        glm::vec3 specular = (F * D * G) / (4 * glm::max(glm::dot(normal, inDir), 0.0f) * glm::max(glm::dot(normal, outDir), 0.0f) + 0.001f); // Adding a small value to avoid division by zero
        // Combine diffuse and specular components
        return diffuse + specular;
    }
    return glm::vec3{0.0f};
  }

  else {
    // Other cases
    // TODO
    return glm::vec3{0.0f};
  }
}

}  // namespace sparks