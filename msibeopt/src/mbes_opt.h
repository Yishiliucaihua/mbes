#pragma once

#include <queue>
#include "./search_interface.h"

using std::queue;
using std::max;
using std::log2;
using std::copy;

namespace core
{
	class mbes_opt : public search_interface
	{
	public:
		mbes_opt(ST _p, ST _q, double _theta, ST _sf, bool _pf, bool _sp, info_reporter *const _ir,
			index_opt *const _idx) : search_interface(_p, _q, _theta, _sf, _pf, _sp, _ir, _idx) {};
		void reduction_1hop(bigraph&);
		void reduction_2hop(bigraph&);
		void solve(const bigraph&);
	private:
		typedef array<pair<ST, ST>, 2> pos_info;
		queue<VIT> tempQ;
		vector<array<ST, 3>> tempC;
		vertex_index_list tempL, tempR, sl, sr, scr;
		unordered_set<VIT> tempUS;
		vs_list tempVS1, tempVS2;
		unordered_map<VIT, ST> l, r;
		vs_list cl, cr;
		ULL solve_imp(const pos_info&, const ST);
		bool intersection_count(const vertex_index_list&, const vs_list&, ST, ST, const ST);
		ST intersection_critical(const vertex_index_list&, const vs_list&, ST, ST, const ST, bool);
		void intersection_ext(const pn_neighbors&, const vs_list&, ST, ST, bool=true);
		void collect(const vs_list&, const ST, const ST);
		void collect_ext(vs_list&, const ST, const ST, vertex_index_list&);
		void update_order(vs_list&, ST, const ST, ST&);
		void restore_order(vs_list&, ST, ST, ST, ST);
		void batch(const vs_list&, const ST, const ST, const vs_list&, const ST, const ST);
		void intersection_lr(const pn_neighbors&, unordered_map<VIT, ST>&);
		void intersection_clcr(const pn_neighbors&, vs_list&, ST, ST);
		void intersection_lr_r(const pn_neighbors&, unordered_map<VIT, ST>&);
		void intersection_clcr_r(const pn_neighbors&, vs_list&, ST, ST);
		void move_clcr_lr(vs_list&, const ST, const ST, unordered_map<VIT, ST>&);
		void update_lr(const unordered_set<VIT>&, unordered_map<VIT, ST>&);
		void update_clcr(const unordered_set<VIT>&, vs_list&, ST, ST);
		void update_lr_r(const vertex_index_list&, unordered_map<VIT, ST>&);
		void update_clcr_r(const vertex_index_list&, vs_list&, ST, ST);
		bool pruning_lr(const unordered_map<VIT, ST>&, const vs_list&, const ST, const ST, const ST, const ST, bool=true);
		void pruning_clcr(vs_list&, const ST, const ST, const vs_list&, const ST, const ST, const ST, const ST);
		bool is_nd(const ST, const vs_list&, ST, ST, const ST);
		bool is_special(const unordered_map<VIT, ST>&, const unordered_map<VIT, ST>&, const vs_list&, const ST, const ST, const vs_list&, const ST, const ST);
		void intersection_special(const pn_neighbors&, const vertex_index_list&, ST, ST, vertex_index_list&);
		void difference_special(const vertex_index_list&, ST, ST, const vertex_index_list&, vertex_index_list&);
		bool pruning_special(const pn_neighbors&, const vertex_index_list&);
		void special(const vertex_index_list&, vertex_index_list&, const vertex_index_list&, vertex_index_list, const ST, const ST);
	};
}
