#pragma once

#include <unordered_map>
#include <algorithm>
#include "./utils.h"

using std::unordered_map;
using std::stable_sort;

namespace core
{
	class index_opt
	{
	public:
		index_opt() {};
		ULL size() const
		{
			return this->s;
		};
		ULL attempts() const
		{
			return this->iters;
		};
		void init(ST);
		void insert(const unordered_map<VIT, ST>&, const unordered_map<VIT, ST>&);
		void insert(const vertex_index_list&, const vertex_index_list&);
		bool insert(const vertex_index_list&);
		void add(const vertex_index_list&);
		void remove_duplicate();
	private:
		ULL s = 0, iters = 0;
		vector<vector<ULL>> data;
		vertex_index_list temp;
		vector<vertex_index_list> vil_data;
		index_opt(const index_opt&) = delete;
		index_opt& operator=(const index_opt&) = delete;
	};
}
