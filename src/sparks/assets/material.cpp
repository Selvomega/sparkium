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
    if (SanityCheck(inDir, outDir, hit_record)) {
      return albedo_color * glm::vec3{scene->GetTextures()[albedo_texture_id].Sample(hit_record.tex_coord)} * float(1/M_PI);
    }
    return glm::vec3{0.0f};
  }
  else if (material_type==MATERIAL_TYPE_EMISSION) {
    // If the material is emissive.
    // Dealt like Lambertian material temporarily but without adding albedo color. 
    // `inDir` is not used here, so it is still fine even if invalid `hit_record.prev_direction` is passed in when explicitly sampling the light. 
    // TODO : Emission BRDF
    std::cout << "EMISSION" << std::endl;
    std::cout << "emission" << emission[0]<< "  " << emission[1]<< "  " << emission[2]<< "  " << std::endl;
    std::cout << "emission_strength" << emission_strength << std::endl;
    return emission * emission_strength;
  }
  else if (material_type == MATERIAL_TYPE_METALLIC) {
    // Cook-Torrance BRDF for metallic surfaces
    std::cout << "METALLIC" << std::endl;
    if (SanityCheck(inDir, outDir, hit_record)) {
      std::cout << "METALLIC aaa" << std::endl;
      glm::vec3 normal = hit_record.normal;
      glm::vec3 halfway = glm::normalize(-inDir + outDir);

      // Specular reflection (Fresnel-Schlick approximation)
      float cosTheta = glm::clamp(glm::dot(normal, halfway), 0.0f, 1.0f);
      // float cosTheta = glm::clamp(-glm::dot(normal, inDir), 0.0f, 1.0f);
      // glm::vec3 F0 = glm::vec3(0.04f); // Assuming non-metallic surface
      glm::vec3 F0 = albedo_color; 
      std::cout << "F0"<< "  " << F0[0]<< "  " << F0[1]<< "  " << F0[2]<< "  " << std::endl;
      std::cout << "cosTheta"<< "  " << cosTheta << std::endl;
      glm::vec3 fresnel = F0 + (1.0f - F0) * pow(1.0f - cosTheta, 5.0f);
      // fresnel become (1,1,1), why?
      std::cout << "fresnel"<< "  " << fresnel[0]<< "  " << fresnel[1]<< "  " << fresnel[2] << std::endl;
      // Microfacet Distribution (Trowbridge-Reitz GGX) 
      float roughness = 0.1; // Roughness value of the surface
      float alpha = roughness * roughness;
      float NdotH = glm::max(glm::dot(hit_record.normal, halfway), 0.0f);
      float D = alpha * alpha / (M_PI * pow(NdotH * NdotH * (alpha * alpha - 1.0f) + 1.0f, 2));
      // float D = 1; // GGX Distribution function

      // Geometric Attenuation (Schlick-GGX approximation)
      float k = (roughness + 1.0f) * (roughness + 1.0f) / 8.0f;
      float NdotL = glm::max(glm::dot(hit_record.normal, inDir), 0.0f);
      float NdotV = glm::max(glm::dot(hit_record.normal, outDir), 0.0f);
      NdotL = glm::max(NdotL, -glm::dot(hit_record.normal, inDir));
      NdotV = glm::max(NdotV, -glm::dot(hit_record.normal, outDir));

      float G1 = NdotL / (NdotL * (1.0f - k) + k);
      float G2 = NdotV / (NdotV * (1.0f - k) + k);
      float G = G1 * G2;
      std::cout << G1 << G2 << G << std::endl;

      std::cout << "D"<< "  " << D << std::endl;
      std::cout << "G"<< "  " << G << std::endl;
      // Final BRDF
      glm::vec3 numerator = fresnel * D * G;
      std::cout << "numerator"<< "  " << numerator[0]<< "  "<< numerator[1]<< "  " << numerator[2]<< "  " << std::endl;
      float denominator = 4 * glm::max(glm::dot(normal, inDir), 0.0f) * glm::max(glm::dot(normal, outDir), 0.0f) + 0.001f; // Adding a small value to avoid division by zero

      glm::vec3 specular = numerator / denominator;

      // Combine diffuse and specular
      std::cout << "specular"<< "  "<<specular[0] << "  "<<specular[1]<< "  " << specular[2]<< "  " << std::endl;
      float max_specular = glm::max(specular[0], glm::max(specular[1], specular[2]));
      std::cout << "max_specular" << max_specular << std::endl;
      // specular = specular / max_specular;
      specular = specular * 2.0f;
      return specular; 
    }
    return glm::vec3{0.0f};
  }
  else if (material_type == MATERIAL_TYPE_DIELECTRIC_GLOSSY) {
    // If the material is a glossy dielectric (non-metallic with specular highlights).
    if (SanityCheck(inDir, outDir, hit_record)) {
        std::cout << "GLOSSY" << std::endl;
        // Albedo color can be either a constant or obtained from a texture
        glm::vec3 albedo = albedo_color; // Assume albedo_color is a glm::vec3 representing the base color
        if (albedo_texture_id >= 0) {
            // If there is a texture, use it to modulate the albedo
            albedo *= glm::vec3{scene->GetTextures()[albedo_texture_id].Sample(hit_record.tex_coord)};
        }

        // Lambertian (Diffuse) component
        glm::vec3 diffuse = albedo * float(1.0 / M_PI);

        // Specular (Cook-Torrance) component
        glm::vec3 halfway = glm::normalize(inDir + outDir);
        float roughness = 0.2; // Roughness value of the surface
        float alpha = roughness * roughness;
        float NdotH = glm::max(glm::dot(hit_record.normal, halfway), 0.0f);
        float D = alpha * alpha / (M_PI * pow(NdotH * NdotH * (alpha * alpha - 1.0f) + 1.0f, 2));

        float k = (roughness + 1.0f) * (roughness + 1.0f) / 8.0f;
        float NdotL = glm::max(glm::dot(hit_record.normal, inDir), 0.0f);
        float NdotV = glm::max(glm::dot(hit_record.normal, outDir), 0.0f);
        NdotV = glm::max(NdotV, -glm::dot(hit_record.normal, outDir));
        float G1 = NdotL / (NdotL * (1.0f - k) + k);
        float G2 = NdotV / (NdotV * (1.0f - k) + k);
        float G = G1 * G2;
        std::cout << 'G' << G1 << G2 << G << std::endl;
        glm::vec3 F0 = glm::vec3(0.04f); // Fresnel at normal incidence for non-metals; can be overridden based on material properties
        float cosTheta = glm::clamp(glm::dot(halfway, outDir), 0.0f, 1.0f);
        glm::vec3 F = F0 + (glm::vec3(1.0f) - F0) * pow(1.0f - cosTheta, 5.0f);

        std::cout << 'F' << F[0] << F[1] << F[2] << std::endl;
        std::cout << 'D' << D << std::endl;
        std::cout << 'G' << G << std::endl;

        glm::vec3 specular = (F * D * G) / (4 * glm::max(glm::dot(hit_record.normal, inDir), 0.0f) * glm::max(glm::dot(hit_record.normal, outDir), 0.0f) + 0.001f); // Adding a small value to avoid division by zero
        std::cout << "specular" << specular[0] << specular[1] << specular[2] << std::endl;
        std::cout << "diffuse" << diffuse[0] << diffuse[1] << diffuse[2] << std::endl;
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