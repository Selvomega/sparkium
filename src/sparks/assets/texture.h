#pragma once
#include "glm/glm.hpp"
#include "string"
#include "vector"

namespace sparks {

// Different types of sampling strategies for the texture. 
enum SampleType { SAMPLE_TYPE_LINEAR = 0, SAMPLE_TYPE_NEAREST = 1 };

class Texture {
 public:
  Texture(uint32_t width = 1,
          uint32_t height = 1,
          const glm::vec4 &color = glm::vec4{1.0f},
          SampleType sample_type = SAMPLE_TYPE_LINEAR);
  Texture(uint32_t width,
          uint32_t height,
          const glm::vec4 *color_buffer,
          SampleType sample_type);
  void Resize(uint32_t width, uint32_t height);
  static bool Load(const std::string &file_path, Texture &texture);
  void Store(const std::string &file_path);
  void SetSampleType(SampleType sample_type);
  [[nodiscard]] SampleType GetSampleType() const;
  glm::vec4 &operator()(int x, int y);
  const glm::vec4 &operator()(int x, int y) const;
  [[nodiscard]] glm::vec4 Sample(glm::vec2 tex_coord) const;
  [[nodiscard]] uint32_t GetWidth() const;
  [[nodiscard]] uint32_t GetHeight() const;
  glm::vec4 *GetBuffer();
  [[nodiscard]] const glm::vec4 *GetBuffer() const;

 private:
  uint32_t width_{};
  uint32_t height_{};
  std::vector<glm::vec4> buffer_; // The color of pixels in the texture.
  SampleType sample_type_{SAMPLE_TYPE_LINEAR};
};
}  // namespace sparks
