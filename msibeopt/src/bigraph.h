#pragma once

#include <fstream>
#include <algorithm>
#include <unordered_map>
#include "./utils.h"

using std::ifstream;
using std::ofstream;
using std::stable_sort;
using std::unordered_map;

namespace core
{
	class bigraph
	{
	public:
		enum class sort_method
		{
			none = 0,
			degree_o,
			degree_positive_o,
			degeneracy_o,
			degeneracy_positive_o,
			method_num
		};
		bigraph() {};
		void load(string);
		void upper(vertex_index_list &_vil) const
		{
			_vil = this->uvil;
		};
		void lower(vertex_index_list &_vil) const
		{
			_vil = this->lvil;
		};
		bool is_swapped() const
		{
			return this->vs1s < this->vs2s;
		};
		ST size() const
		{
			return this->va;
		};
		const vertex_set& data() const
		{
			return this->vs;
		};
		void print_info() const;
		void reorder(sort_method);
		void compact();
		void write(string) const;
	private:
		friend class mbes_opt;
		vertex_set vs;
		vertex_index_list uvil, lvil;
		vector<array<unordered_set<VIT>, 2>> cache;
		ST va = 0, ea = 0, pea = 0, vs1s = 0, vs2s = 0;
		void sort_p(sort_method, vs_list&);
	};
}
