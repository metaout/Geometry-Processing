#pragma once
#include <Eigen/Eigen>
#include "halfedge.hpp"
#include "union_find.hpp"

namespace gp {
	constexpr double m_PI = 3.14159265359;

	template <typename T>
	double cot(const T edge) {
		auto e1 = edge->vert->pos.eig - edge->next->vert->pos.eig;
		auto e2 = edge->prev->vert->pos.eig - edge->next->vert->pos.eig;

		double dot = e1.dot(e2);
		double cross_norm = (e1.cross(e2)).norm();

		return dot / cross_norm;
	}

	Eigen::SparseMatrix<double> areaSparseMat(const heds::Mesh& mesh, bool mean = true) {
		using Triplet = Eigen::Triplet<double>;
		Eigen::SparseMatrix<double> mat(mesh.getNumVertices(), mesh.getNumVertices());
		std::vector<Triplet> triplets;

		Eigen::VectorXd count = Eigen::VectorXd::Zero(mesh.getNumVertices());
		for (auto& f : mesh.faces) {
			auto a = f->area();
			auto v0 = f->edge->vert;
			auto v1 = f->edge->next->vert;
			auto v2 = f->edge->next->next->vert;

			triplets.push_back(Triplet(v0->id, v0->id, a));
			triplets.push_back(Triplet(v1->id, v1->id, a));
			triplets.push_back(Triplet(v2->id, v2->id, a));

			count(v0->id) += 1.0;
			count(v1->id) += 1.0;
			count(v2->id) += 1.0;
		}

		mat.setFromTriplets(triplets.begin(), triplets.end());
		if (mean) mat.diagonal().array() /= count.array();
		return mat;
	}

	Eigen::VectorXd areaVectorPerFace(const heds::Mesh& mesh, float div = 3.0) {
		Eigen::VectorXd areas(mesh.getNumFaces());
		for (auto& face : mesh.faces) areas[face->id] = face->area() / div;
		return areas;
	}

	Eigen::VectorXd areaVectorPerVert(const heds::Mesh& mesh, float div = 3.0) {
		Eigen::VectorXd areas = Eigen::VectorXd::Zero(mesh.getNumVertices());
		for (auto& face : mesh.faces) {
			auto a = face->area();
			auto v0 = face->edge->vert;
			auto v1 = face->edge->next->vert;
			auto v2 = face->edge->next->next->vert;
			areas[v0->id] += a;
			areas[v1->id] += a;
			areas[v2->id] += a;
		}
		areas.array() /= div;
		return areas;
	}

	Eigen::SparseMatrix<double> cotanLaplacianSparseMat(const heds::Mesh& mesh,
																											const bool normalize = true) {
		using Triplet = Eigen::Triplet<double>;
		Eigen::SparseMatrix<double> mat(mesh.getNumVertices(), mesh.getNumVertices());
		std::vector<Triplet> triplets;
		
		Eigen::VectorXd sum_weights = Eigen::VectorXd::Zero(mesh.verts.size());
		for (auto& face : mesh.faces) {
			auto v0 = face->edge->vert;
			auto v1 = face->edge->next->vert;
			auto v2 = face->edge->next->next->vert;

			auto w0 = 0.5 * cot(mesh.edges.at({ v1, v2 }));
			auto w1 = 0.5 * cot(mesh.edges.at({ v2, v0 }));
			auto w2 = 0.5 * cot(mesh.edges.at({ v0, v1 }));

			triplets.push_back(Triplet(v0->id, v0->id, w2 + w1));
			triplets.push_back(Triplet(v1->id, v1->id, w2 + w0));
			triplets.push_back(Triplet(v2->id, v2->id, w0 + w1));

			triplets.push_back(Triplet(v0->id, v1->id, -w2));
			triplets.push_back(Triplet(v1->id, v0->id, -w2));
			triplets.push_back(Triplet(v1->id, v2->id, -w0));
			triplets.push_back(Triplet(v2->id, v1->id, -w0));
			triplets.push_back(Triplet(v2->id, v0->id, -w1));
			triplets.push_back(Triplet(v0->id, v2->id, -w1));

			sum_weights[v0->id] += (w2 + w1);
			sum_weights[v1->id] += (w2 + w0);
			sum_weights[v2->id] += (w1 + w0);
		}

		mat.setFromTriplets(triplets.begin(), triplets.end());
		if (normalize) {
			return  areaVectorPerVert(mesh).asDiagonal().inverse() * mat;
		}
		return mat;
	}

	Eigen::VectorXd gaussianCurvature(heds::Mesh& mesh) {
		Eigen::VectorXd curvature = Eigen::VectorXd::Zero(mesh.getNumVertices());
		curvature.fill(2.0 * m_PI);

		for (auto& face : mesh.faces) {
			auto v1 = face->edge->vert;
			auto v2 = face->edge->next->vert;
			auto v3 = face->edge->next->next->vert;

			auto l12 = v2->pos.eig - v1->pos.eig;
			auto l13 = v3->pos.eig - v1->pos.eig;
			auto l21 = v1->pos.eig - v2->pos.eig;
			auto l23 = v3->pos.eig - v2->pos.eig;
			auto l31 = v1->pos.eig - v3->pos.eig;
			auto l32 = v2->pos.eig - v3->pos.eig;

			auto theta1 = acos(l12.dot(l13) / (l12.norm() * l13.norm()));
			auto theta2 = acos(l21.dot(l23) / (l21.norm() * l23.norm()));
			auto theta3 = acos(l31.dot(l32) / (l31.norm() * l32.norm()));

			curvature[v1->id] -= theta1;
			curvature[v2->id] -= theta2;
			curvature[v3->id] -= theta3;
		}
		
		auto area = areaVectorPerVert(mesh);
		curvature.array() /= area.array();
		return curvature;
	}

	Eigen::VectorXd meanCurvature(heds::Mesh& mesh, Eigen::SparseMatrix<double>& lbo) {
		Eigen::VectorXd curvature = Eigen::VectorXd::Zero(mesh.getNumVertices());
		Eigen::MatrixXd verts = Eigen::MatrixXd(mesh.getNumVertices(), 3);

		for (size_t i = 0; i < mesh.getNumVertices(); i++) verts.row(i) = mesh.verts[i]->pos.eig;
		Eigen::MatrixXd H = lbo * verts;
		for (size_t i = 0; i < mesh.getNumVertices(); i++) curvature(i) = 0.5 * H.row(i).norm();
		
		return curvature;
	}

	Eigen::VectorXd heatDiffusionSmoothing(Eigen::SparseMatrix<double>& lbo, Eigen::VectorXd heat, const double time) {
		Eigen::SparseMatrix<double> eye(lbo.rows(), lbo.rows());
		eye.setIdentity();
		Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
		solver.compute(eye + time * lbo);

		Eigen::VectorXd smoothed_heat = solver.solve(heat);
		return smoothed_heat;
	}

	template <typename T1, typename T2>
	auto connectedMesh(const T1& pos, const T2& faces) {
		UnionFind uf(pos.size());
		for (auto& face : faces) {
			for (size_t i = 1; i < face.size(); i++) uf.unite(face[0], face[i]);
		}
		
		std::unordered_map<size_t, std::vector<size_t>> s;
		for (size_t i = 0; i < pos.size(); i++) {
			s[uf.root(i)].push_back(i);
		}

		return s;
	}
};