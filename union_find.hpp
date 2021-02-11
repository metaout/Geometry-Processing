#pragma once
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>
#include <unordered_set>

namespace gp {
	class UnionFind {
		size_t numNode;
		std::vector<size_t> parents;

	public:
		UnionFind(size_t n) : numNode(n), parents(n) {
			std::iota(parents.begin(), parents.end(), 0);
		}

		size_t root(size_t node) {
			if (node == parents[node]) {
				return node;
			} else {
				return parents[node] = root(parents[node]);
			}
		}

		bool unite(size_t lnode, size_t rnode) {
			size_t lr = root(lnode), rr = root(rnode);
			if (lr == rr) return false;

			parents[lr] = rr;
			return true;
		}

		bool same(size_t lnode, size_t rnode) {
			return root(lnode) == root(rnode);
		}
	};
}