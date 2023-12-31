#pragma once
#include "sparks/assets/aabb.h"
#include "sparks/assets/model.h"
#include "sparks/assets/util.h"
#include "sparks/assets/vertex.h"
#include "sparks/util/mikktspace.h"
#include "vector"


namespace sparks {
class Mesh : public Model {
 public:
  Mesh() = default;
  Mesh(const Mesh &mesh);
  Mesh(const std::vector<Vertex> &vertices,
       const std::vector<uint32_t> &indices);
  explicit Mesh(const tinyxml2::XMLElement *element);
  ~Mesh() override = default;
  [[nodiscard]] float TraceRay(const glm::vec3 &origin,
                               const glm::vec3 &direction,
                               float t_min,
                               HitRecord *hit_record) const override;
  const char *GetDefaultEntityName() override;
  [[nodiscard]] AxisAlignedBoundingBox GetAABB(
      const glm::mat4 &transform) const override;
  [[nodiscard]] std::vector<Vertex> GetVertices() const override;
  [[nodiscard]] std::vector<uint32_t> GetIndices() const override;
  static Mesh Cube(const glm::vec3 &center, const glm::vec3 &size);
  static Mesh Sphere(const glm::vec3 &center = glm::vec3{0.0f},
                     float radius = 1.0f);
  static bool LoadObjFile(const std::string &obj_file_path, Mesh &mesh);
  void WriteObjFile(const std::string &file_path) const;
  void MergeVertices();
  glm::vec3 Random() const override;

 //protected:
  std::vector<Vertex> vertices_;
  std::vector<uint32_t> indices_;
  AxisAlignedBoundingBox boundingbox_;
  void InitBoundingBox();

  static int getNumFaces(const SMikkTSpaceContext *pContext);
  static int getNumVerticesOfFace(const SMikkTSpaceContext *pContext, const int iFace);
  static void getPosition(const SMikkTSpaceContext *pContext, float fvPosOut[], const int iFace, const int iVert);
  static void getNormal(const SMikkTSpaceContext *pContext, float fvNormOut[], const int iFace, const int iVert);
  static void getTexCoord(const SMikkTSpaceContext *pContext, float fvTexcOut[], const int iFace, const int iVert);
  static void setTSpaceBasic(const SMikkTSpaceContext *pContext, const float fvTangent[], const float fSign, const int iFace, const int iVert);
};
}  // namespace sparks
