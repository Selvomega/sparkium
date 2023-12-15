#include "sparks/assets/accelerated_mesh.h"
#include <cassert>
#include <iterator>
#include <vector>

#include "algorithm"
#include "sparks/assets/aabb.h"
#include "sparks/assets/mesh.h"
#include "sparks/assets/vertex.h"

namespace sparks {
BVHTree::BVHTree() {
	root = nullptr;
}
BVHTree::BVHTree(std::vector<AxisAlignedBoundingBox> &b) {
	std::vector<std::pair<const AxisAlignedBoundingBox *,int>> objs;
	buf = std::move(b); // 'b' cannot be used anymore. 
	buf.shrink_to_fit(); // To save space.
	objs.reserve(buf.size());
	const AxisAlignedBoundingBox *ptr = &buf[0];
	for (size_t i = 0; i < buf.size(); i++) {
		objs.emplace_back(std::make_pair(&ptr[i],i));
	}
	root.reset(Build(objs));
}

BVHTree::BVHNode *
BVHTree::Build(std::vector<std::pair<const AxisAlignedBoundingBox *,int>> &objs, int depth) {
	/*
	Build the BVH tree. 
	The tree guarantees that there is no node residing in different branches. 
	What is `depth` used for?
	*/
	BVHTree::BVHNode *root = new BVHTree::BVHNode();
	auto& aabb = root->aabb;
	aabb = *(objs[0].first);
	for (auto ptr:objs) {
    // Compute the total bounding box. 
    aabb = aabb | *(ptr.first);
  }
	if (objs.size() < 3) {
    // 3 is the magic number, each aabb can contain at most 4 objects. 
		root->left = nullptr;
		root->right = nullptr;
		root->objs = objs;
		return root;
	} 
  // Assume the longest distance is along the x-axis
	int dim = 2;
	auto x_distance = aabb.x_high-aabb.x_low;
	auto y_distance = aabb.y_high-aabb.y_low;
  auto z_distance = aabb.z_high-aabb.z_low;
	if (x_distance > y_distance && x_distance > z_distance) {
    // If the longest distance is along the x-axis
		dim = 0;
  }
	else if (y_distance > z_distance) {
    // If the longest distance is along the y-axis
		dim = 1;
  }

  std::sort(objs.begin(), objs.end(), [dim](auto a, auto b) {
    return ((a.first)->centroid(dim)) < ((b.first)->centroid(dim));
	});
  
  // Chop along the longest direction. 
	auto b = objs.begin();
	auto m = b + (objs.size() / 2);
	auto e = objs.end();

	auto left = std::vector<std::pair<const AxisAlignedBoundingBox *,int>>(b, m);
	auto right = std::vector<std::pair<const AxisAlignedBoundingBox *,int>>(m, e);
	root->left.reset(Build(left, depth + 1));
	root->right.reset(Build(right, depth + 1));
	return root;
}

void 
BVHTree::Intersect_r(BVHNode *root,
                     const glm::vec3 &origin,
                     const glm::vec3 &direction,
                     float t_min,
                     std::vector<int> &result) const {
	if (root==nullptr) {
		return;
	}
	if (!root->aabb.IsIntersect(origin, direction, t_min, 1e5)) {
		// If the ray has no intersection with the AABB. 
		return;
	}
	// If the ray hits the AABB. 
	if (root->left.get() == nullptr) {
		// When the left child of the node is empty. 
		// Then the right child must also be empty. 
		// Which means that this is a leaf, with at most 5 objects in it. 
		assert(root->right.get() == nullptr);
		// const AxisAlignedBoundingBox *start = &buf[0]; // Used to compute the offset. 
		for (auto ptr:root->objs) {
			// Only when we have reached the leaf will the objects be added to the `result` vector. 
			// So no objects will be counted twice. 
			result.emplace_back(ptr.second);
		}
//std::cout<<"Hit size: "<<result.size()<<std::endl;
		return;
	}
	Intersect_r(root->left.get(), origin, direction, t_min, result);
	Intersect_r(root->right.get(), origin, direction, t_min, result);
}

bool
BVHTree::Intersect(const glm::vec3 &origin,
                   const glm::vec3 &direction,
                   float t_min,
                   std::vector<int> &result) const {
	result.clear();
	Intersect_r(root.get(), origin, direction, t_min, result);
	return result.size(); // Whether there are intersections. 
}

AxisAlignedBoundingBox 
BVHTree::AccessBuf(int index) const{
	return buf[index];
}

AcceleratedMesh::AcceleratedMesh(const Mesh &mesh) : Mesh(mesh) {
	BuildAccelerationStructure();
}

AcceleratedMesh::AcceleratedMesh(const std::vector<Vertex> &vertices,
                                 const std::vector<uint32_t> &indices)
    : Mesh(vertices, indices) {
		BuildAccelerationStructure();
}

float AcceleratedMesh::TraceRay(const glm::vec3 &origin,
                                const glm::vec3 &direction,
                                float t_min,
                                HitRecord *hit_record) const {
	float ret = -1;
	std::vector<int> result={};
	bool not_miss = tree.Intersect(origin, direction, t_min, result);
	if (!not_miss) {
		return ret;
	}
	assert(result.size()>0);
	for (int index : result) {
		AxisAlignedBoundingBox aabb = tree.AccessBuf(index);
		for (int i=0; i<aabb.vertices_.size(); i+=3) {
			int j = i + 1, k = i + 2;
			const auto &v0 = aabb.vertices_[i];
			const auto &v1 = aabb.vertices_[j];
			const auto &v2 = aabb.vertices_[k];
			glm::mat3 A = glm::mat3(v1.position - v0.position,
                              v2.position - v0.position, -direction);
			if (std::abs(glm::determinant(A)) < 1e-9f) {
				// The direction is parallel to the plane of the triangle. 
				continue;
			}
			A = glm::inverse(A);
			auto uvt = A * (origin - v0.position);
			// `uvt.x` and `uvt.y` determine the intersection point of the ray with the plane of the triangle. 
			auto &t = uvt.z;
			if (t < t_min || (ret > 0.0f && t > ret)) {
				// If the distance is too close or there has already been closer intersections. 
				continue;
			}
			auto &u = uvt.x;
			auto &v = uvt.y;
			auto w = 1.0f - u - v;
			auto position = origin + t * direction; // The intersection position. 
			if (u >= 0.0f && v >= 0.0f && u + v <= 1.0f) {
				// The intersection point lies in the triangle. 
				ret = t;
				if (hit_record) {
					auto geometry_normal = glm::normalize(
							glm::cross(v2.position - v0.position, v1.position - v0.position));
					// The geometry normal is the normal of model. 
					if (glm::dot(geometry_normal, direction) < 0.0f) {
						hit_record->position = position;
						hit_record->geometry_normal = geometry_normal; // To make sure that the geometry normal points to the direction where the ray comes.
						hit_record->normal = v0.normal * w + v1.normal * u + v2.normal * v; // The normal is the processed normal. 
						hit_record->tangent =
								v0.tangent * w + v1.tangent * u + v2.tangent * v; // What is the weighted tagent used for?
						hit_record->tex_coord =
								v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
						hit_record->front_face = true;
					} else {
						hit_record->position = position;
						hit_record->geometry_normal = -geometry_normal;
						hit_record->normal = -(v0.normal * w + v1.normal * u + v2.normal * v);
						hit_record->tangent =
								-(v0.tangent * w + v1.tangent * u + v2.tangent * v);
						hit_record->tex_coord =
								v0.tex_coord * w + v1.tex_coord * u + v2.tex_coord * v;
						hit_record->front_face = false;
					}
				}
			}
		}
	}
	return ret;
}

void AcceleratedMesh::BuildAccelerationStructure() {
	std::vector<AxisAlignedBoundingBox> vec;
	for (int i=0; i<indices_.size(); i+=3) {
		std::vector<Vertex> temp_vec = {};
		int j = i + 1, k = i + 2;
		temp_vec.push_back(vertices_[indices_[i]]);
		temp_vec.push_back(vertices_[indices_[j]]);
		temp_vec.push_back(vertices_[indices_[k]]);
    vec.push_back(AxisAlignedBoundingBox(temp_vec));
	}
	tree = BVHTree(vec);
}
// TODO
}  // namespace sparks
