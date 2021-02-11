#pragma once
#include <algorithm>
#include <utility>
#include <vector>
#include <memory>

namespace gp {
	template<typename T>
	class kdTree {
		std::vector<T> points;
		int depth;
		int dim;
		int median;
		std::shared_ptr<kdTree<T>> left, right;

	public:
		kdTree(std::vector<T> points, int depth = 0, int dim = 3) : points(points), depth(depth), dim(dim), median(0),
			left(nullptr), right(nullptr) {}

		void subdiv() {
			if (left == nullptr) {
				int axis = this->depth % this->dim;

				std::sort(points.begin(), points.end(), [&](auto l, auto r) { return l[axis] < r[axis]; });
				median = points.size() / 2;
				left = std::make_shared<kdTree>(std::vector<T>(points.begin(), points.begin() + median), depth + 1, dim);
				right = std::make_shared<kdTree>(std::vector<T>(points.begin() + median, points.end()), depth + 1, dim);

				return;
			}
			left->subdiv();
			right->subdiv();
		}

		std::vector<std::vector<T>> getChildPoints(std::shared_ptr<kdTree<T>>& tree) {
			if (tree->left == nullptr) return { tree->points };

			std::vector<std::vector<T>> p;
			std::vector<std::vector<T>> lp = getChildPoints(tree->left);
			std::vector<std::vector<T>> rp = getChildPoints(tree->right);

			p.insert(p.end(), lp.begin(), lp.end());
			p.insert(p.end(), rp.begin(), rp.end());

			return p;
		}

		std::vector<std::vector<T>> getPoints() {
			if (left == nullptr) return { points };

			std::vector<std::vector<T>> p;
			std::vector<std::vector<T>> lp = getChildPoints(left);
			std::vector<std::vector<T>> rp = getChildPoints(right);

			p.insert(p.end(), lp.begin(), lp.end());
			p.insert(p.end(), rp.begin(), rp.end());

			return p;
		}

		std::vector<T> getNearestNeighbors(T p) {
			if (left == nullptr) return points;

			int axis = this->depth % this->dim;
			auto tree = points[this->median][axis] > p[axis] ? left : right;
			while (tree->left != nullptr) {
				axis = tree->depth % this->dim;
				if (tree->points[tree->median][axis] > p[axis]) tree = tree->left;
				else tree = tree->right;
			}

			return tree->points;
		}
	};
}