#include "sparks/renderer/path_tracer.h"
#include <iostream>
#include <queue>

#include "glm/fwd.hpp"
#include "glm/geometric.hpp"
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
  glm::vec3 prev_throughput{1.0f};
  glm::vec3 radiance{0.0f};
  HitRecord hit_record;
  const int max_bounce = render_settings_->num_bounces;
  std::mt19937 rd(sample ^ x ^ y); // Pseudorandom generator. 
  glm::vec3 accumulated{1.0f};
  for (int j = 0; j < max_bounce; j++) {
    for (int i=0; i<scene_->GetEntityCount(); i++) {
      // Sample the light explicitly. 
      // Is it too slow? 
      if (hit_record.hit_entity_id==-1) {
        break;
      }
      if (scene_->GetEntity(i).GetMaterial().material_type==MATERIAL_TYPE_EMISSION && i!=hit_record.hit_entity_id) {
        // If the material is emissive. 
        const Mesh* mesh = dynamic_cast<const Mesh*>(scene_->GetEntity(i).GetModel());
        const glm::mat4 transform = scene_->GetEntity(i).GetTransformMatrix();
        for (int first=0; first<mesh->indices_.size(); first+=3) {
          // Enumerate the triangles of the mesh. 
          int second = first + 1, third = first + 2;
          const auto &v0 = mesh->vertices_[mesh->indices_[first]];
          const auto &v1 = mesh->vertices_[mesh->indices_[second]];
          const auto &v2 = mesh->vertices_[mesh->indices_[third]];
          std::random_device nrd;
          std::mt19937 gen(nrd());
          std::uniform_real_distribution<> dis(0.0, 1.0);
          float u_rand = dis(gen);
          float v_rand = dis(gen);
          float w_rand = dis(gen);
          float sum = u_rand+v_rand+w_rand;
          u_rand = u_rand / sum;
          v_rand = v_rand / sum;
          w_rand = w_rand / sum;
          auto mesh_rand_position = u_rand*v0.position + v_rand*v1.position + w_rand*v2.position;
          glm::vec3 global_rand_position = transform*glm::vec4{mesh_rand_position,1.0};
          glm::vec3 raw_direction = global_rand_position-origin;
          if (glm::length(raw_direction)<1e-3 || (hit_record.hit_entity_id!=-1 && (glm::dot(raw_direction,hit_record.geometry_normal)<0))) {
            // If the two points are too close. 
            // Or the mesh piece is on the other side of the normal.
            // You can believe that the geometry normal points to the outside of the object. 
            continue; // Next triangle.
          }
          glm::vec3 global_direction = glm::normalize(raw_direction);
          HitRecord temp_hit_record;
          auto t = scene_->TraceRay(origin, global_direction, 1e-3f, 1e4f, &temp_hit_record);
          if (t<=0.0 || std::abs(t-glm::length(raw_direction))>1e-2) {
            // If did not hit the light source. 
            continue; // Next triangle.
          }
          auto material = scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
          auto light_material = scene_->GetEntity(i).GetMaterial();
          auto cosine1 = std::abs(glm::dot(global_direction, temp_hit_record.normal)); 
          auto cosine2 = std::abs(glm::dot(global_direction, hit_record.normal)); 
          radiance += prev_throughput*light_material.emission*light_material.emission_strength*material.BRDF(hit_record.prev_direction, direction, hit_record, scene_) * cosine1 * cosine2 * glm::length(glm::cross(v1.position-v0.position, v2.position-v0.position)) / (2*glm::length(raw_direction)*glm::length(raw_direction));
        }
      }
    }
    // `hit_record` is modified at this step. 
    auto t = scene_->TraceRay(origin, direction, 1e-3f, 1e4f, &hit_record);
    if (t > 0.0f) {
      // Hit an object. 
      auto &material =
          scene_->GetEntity(hit_record.hit_entity_id).GetMaterial();
      if (material.material_type==MATERIAL_TYPE_EMISSION) {
        // If the material is emissive. 
        if (j==0) {
          // If the first radiance
          radiance += throughput * material.emission * material.emission_strength;
        }
        break; // Break since we have explicitly sampled the light source. 
      }
      else if (material.material_type==MATERIAL_TYPE_SPECULAR) {
        // If the material is totally specular. 
        origin = hit_record.position;
        hit_record.prev_direction = direction;
        // The following line is true no matter what is the direction of the normal
        direction = direction - 2*glm::dot(direction,hit_record.normal)*hit_record.normal;
        if (!SanityCheck(hit_record.prev_direction, direction, hit_record)) {
          // If sanity check not passed, 
          // Eliminate the ray? 
          radiance *= 0;
          break;
        }
        // Specular reflection, continue
      }
      else if (material.material_type==MATERIAL_TYPE_LAMBERTIAN) {
        // Lambertian material 
        origin = hit_record.position;
        auto p = material.ImportanceSampling(direction, hit_record);
        auto out_direction = p.first;
        auto pdf = p.second;
        // This line?
        auto cosine = std::abs(glm::dot(out_direction, hit_record.normal)); 
        prev_throughput = throughput;
        throughput *= material.BRDF(direction, out_direction, hit_record, scene_)*cosine/pdf;
        hit_record.prev_direction = direction;
        direction = out_direction;
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
      /*
      if (j > 3) {
        float p = std::max(throughput.x, std::max(throughput.y, throughput.z));
        std::random_device nrd;
        std::mt19937 gen(nrd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        float rand_x = dis(gen);
        if (rand_x > p) {
          break;
        }
        throughput *= 1 / p;
      }
      */
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
