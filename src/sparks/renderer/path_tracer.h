#pragma once
#include "random"
#include "sparks/assets/scene.h"
#include "sparks/renderer/renderer_settings.h"

namespace sparks {
class PathTracer {
 public:
  PathTracer(const RendererSettings *render_settings, const Scene *scene);
  [[nodiscard]] glm::vec3 SampleRay(glm::vec3 origin,
                                    glm::vec3 direction,
                                    int x,
                                    int y,
                                    int sample);

  [[nodiscard]] glm::vec3 SampleRayOld(glm::vec3 origin,
                                    glm::vec3 direction,
                                    int x,
                                    int y,
                                    int sample) const;

  [[nodiscard]] static std::pair<glm::vec3,float> RandomSampling(const glm::vec3 &inDir, const HitRecord &hit_record);
  [[nodiscard]] static std::pair<glm::vec3,float> ImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record);
  [[nodiscard]] static std::pair<glm::vec3,float> CosImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record);
  [[nodiscard]] static std::pair<glm::vec3,float> MultiImportanceSampling(const glm::vec3 &inDir, const HitRecord &hit_record);

  bool GetEntered(int index);
  void SetEntered(int index, bool value); 
  
 private:
  std::vector<bool> entered_;
  const RendererSettings *render_settings_{};
  const Scene *scene_{};
};
}  // namespace sparks
