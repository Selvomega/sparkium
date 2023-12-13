#include "sparks/renderer/path_tracer.h"

#include "sparks/assets/material.h"
#include "sparks/assets/util.h"
#include "sparks/util/util.h"

namespace sparks {
PathTracer::PathTracer(const RendererSettings *render_settings,
                       const Scene *scene) {
  render_settings_ = render_settings;
  scene_ = scene;
}

glm::vec3 PathTracer::SampleRay(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) const {
  // TODO
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  HitRecord hit_record;
  const int max_bounce = render_settings_->num_bounces;
  std::mt19937 rd(sample ^ x ^ y); // Pseudorandom generator. 
  glm::vec3 accumulated{1.0f};
  for (int i = 0; i < max_bounce; i++) {
    // `hit_record` is modified at this step. 
    auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record);
    if (t > 0.0f) {
      // Hit an object. 
      auto &material =
          scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
      if (material.material_type==MATERIAL_TYPE_EMISSION) {
        // If the material is emissive. 
        radiance += throughput * material.emission * material.emission_strength;
        break;
      }
      else if (material.material_type==MATERIAL_TYPE_SPECULAR) {
        // If the material is totally specular. 
        origin = hit_record.position;
        direction = direction - 2*glm::dot(direction,hit_record.normal)*hit_record.normal;
        if (!SanityCheck(direction, hit_record)) {
          // If sanity check not passed, 
          // Eliminate the ray? 
          radiance *= 0;
          break;
        }
        // Specular reflection, continue
        continue; 
      }
      else if (material.material_type==MATERIAL_TYPE_LAMBERTIAN) {
        // Lambertian material 
        origin = hit_record.position;
        auto p = material.ImportanceSampling(direction, hit_record);
        auto out_direction = p.first;
        auto pdf = p.second;
        auto cosine = glm::dot(out_direction, hit_record.normal); 
        throughput *= material.BRDF(direction, out_direction, hit_record, scene_)*cosine/pdf;
        direction = out_direction;
        continue;
      }
      else {
        // The other cases. 
        // TODO
        throughput *=
            material.albedo_color *
            glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(
                hit_record.tex_coord)};
        origin = hit_record.position;
        // Recall that the environment is a huge sphere at "infinity"
        direction = scene_->GetEnvmapLightDirection();
        radiance += throughput * scene_->GetEnvmapMinorColor();
        throughput *=
            std::max(glm::dot(direction, hit_record.normal), 0.0f) * 2.0f;
        if (scene_->TraceRay(origin, direction, 1e-3f, 1e4f, nullptr) < 0.0f) {
          // If the ray to the environment is not blocked, which means that the environment light can directly reach the object. 
          radiance += throughput * scene_->GetEnvmapMajorColor();
        }
        break;
      }
    }
    else {
      // Did not hit any object. 
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(direction)};
      break;
    }
  }
  return radiance;
}

glm::vec3 PathTracer::SampleRayOld(glm::vec3 origin,
                                glm::vec3 direction,
                                int x,
                                int y,
                                int sample) const {
  // TODO
  glm::vec3 throughput{1.0f};
  glm::vec3 radiance{0.0f};
  HitRecord hit_record;
  const int max_bounce = render_settings_->num_bounces;
  std::mt19937 rd(sample ^ x ^ y); // Pseudorandom generator. 
  for (int i = 0; i < max_bounce; i++) {
    auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record);
    if (t > 0.0f) {
      // Get the material of the hit variable. 
      auto &material =
          scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
      if (material.material_type == MATERIAL_TYPE_EMISSION) {
        radiance += throughput * material.emission * material.emission_strength;
        break;
      } else {
        throughput *=
            material.albedo_color *
            glm::vec3{scene_->GetTextures()[material.albedo_texture_id].Sample(
                hit_record.tex_coord)};
        origin = hit_record.position;
        // Recall that the environment is a huge sphere at "infinity"
        direction = scene_->GetEnvmapLightDirection();
        radiance += throughput * scene_->GetEnvmapMinorColor();
        throughput *=
            std::max(glm::dot(direction, hit_record.normal), 0.0f) * 2.0f;
        if (scene_->TraceRay(origin, direction, 1e-3f, 1e4f, nullptr) < 0.0f) {
          // If the ray to the environment is not blocked, which means that the environment light can directly reach the object. 
          radiance += throughput * scene_->GetEnvmapMajorColor();
        }
        break;
      }
    } else {
      radiance += throughput * glm::vec3{scene_->SampleEnvmap(direction)};
      break;
    }
  }
  return radiance;
}
}  // namespace sparks
