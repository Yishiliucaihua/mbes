#pragma once

#include <queue>
#include <list>
#include "./search_interface.h"

using std::queue;
using std::list;
using std::max;

namespace core
{
	class mbes : public search_interface
	{
	public:
		mbes(ST _p, ST _q, double _theta, ST _sf, bool _pf, bool _sp, info_reporter *const _ir,
			index_opt *const _idx) : search_interface(_p, _q, _theta, _sf, _pf, _sp, _ir, _idx) {};
		void reduction_1hop(bigraph&);
		void reduction_2hop(bigraph&);
		void solve(const bigraph&);
	private:
		vertex_index_list tempV, txr;
		queue<VIT> tempQ;
		unordered_set<VIT> tempUS, tempS1, tempS2, tl, tr, tcr;
		vector<array<ST, 3>> tempC;
		vs_map tempVM1, tempVM2;
		ULL solve_imp(const vs_map&, const vs_map&, const vs_map&, const vs_map&);
		void intersection_lr(const pn_neighbors&, vs_map&);
		void intersection_clcr(const pn_neighbors&, const vs_map&, vs_map&);
		void update_lr(const unordered_set<VIT>&, vs_map&);
		void update_clcr(const unordered_set<VIT>&, const vs_map&, vs_map&);
		bool pruning_lr(const vs_map&, const ST, const vs_map&, const ST);
		void pruning_clcr(const vs_map&, const ST, const vs_map&, const ST);
		bool is_nd(const ST, const vs_map&, const ST);
		bool is_special(const vs_map&, const vs_map&, const vs_map&, const vs_map&);
		void special(const unordered_set<VIT>&, const unordered_set<VIT>&, const unordered_set<VIT>&, vertex_index_list, const ST, const ST);
	};
}
