#pragma once

#include <cmath>
#include "./bigraph.h"
#include "./index.h"

namespace core
{
	class search_interface
	{
	public:
		search_interface(ST _p, ST _q, double _theta, ST _sf, bool _pf, bool _sp, info_reporter *const _ir,
			index_opt *const _idx) : p(_p), q(_q), pdcL((ST)ceil(_theta * this->q)), pdcR((ST)ceil(_theta * this->p)),
			theta(_theta), omt(1.0f - _theta), sf(_sf), pf(_pf), sp(_sp), ir(_ir), idx(_idx) {};
		virtual void reduction_1hop(bigraph&) = 0;
		virtual void reduction_2hop(bigraph&) = 0;
		virtual void solve(const bigraph&) = 0;
		ULL size() const
		{
			return this->idx->size();
		};
	protected:
		ST p, q, pdcL, pdcR;
		double theta, omt;
		ST sf = 0;  // 0: disable swap; 1: enable swap
		bool pf = false;  // false: disable pruning; true: enable pruning
		bool sp = false;  // false: disable special case; true: enable special case
		info_reporter *const ir;
		index_opt *const idx;
		const vertex_set *vs = nullptr;
		double iters = 0.0;
	};
}
