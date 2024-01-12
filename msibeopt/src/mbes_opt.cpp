#include "./mbes_opt.h"

using namespace core;

void mbes_opt::reduction_1hop(bigraph &_bg)
{
	// clear queue
	while (this->tempQ.size())
	{
		this->tempQ.pop();
	}

	// initialize queue
	const VIT va = _bg.va, vs1s = _bg.vs1s;
	auto &vs = _bg.cache;

	auto init = [&vs, this](const VIT _s, const VIT _e, const ST _pdc, const ST _dc)
	{
		for (VIT i = _s; i < _e; ++i)
		{
			if ((ST)vs[i][0].size() < _pdc || (ST)(vs[i][0].size() + vs[i][1].size()) < _dc)
			{
				this->tempQ.emplace(i);
			}
		}
	};

	init(0, vs1s, this->pdcL, this->q);
	init(vs1s, va, this->pdcR, this->p);

	auto update = [&vs, this](const VIT _vi, const ST _pdc, const ST _dc)
	{
		for (const auto pn : vs[_vi][0])
		{
			vs[pn][0].erase(_vi);
			if ((ST)vs[pn][0].size() < _pdc || (ST)(vs[pn][0].size() + vs[pn][1].size()) < _dc)
			{
				this->tempQ.push(pn);
			}
		}
		for (const auto nn : vs[_vi][1])
		{
			vs[nn][1].erase(_vi);
			if ((ST)(vs[nn][0].size() + vs[nn][1].size()) < _dc)
			{
				this->tempQ.push(nn);
			}
		}
	};

	// peeling
	while (this->tempQ.size())
	{
		VIT vi = this->tempQ.front();
		this->tempQ.pop();

		if (vi < vs1s)
		{
			update(vi, this->pdcR, this->p);
		}
		else
		{
			update(vi, this->pdcL, this->q);
		}

		vs[vi][0].clear();
		vs[vi][1].clear();
	}

	_bg.compact();
}

void mbes_opt::reduction_2hop(bigraph &_bg)
{
	// clear queue
	while (this->tempQ.size())
	{
		this->tempQ.pop();
	}

	// initialize queue
	const VIT va = _bg.va, vs1s = _bg.vs1s;
	auto &vs = _bg.cache;
	this->tempC.resize(va);

	auto init = [&vs, this](const VIT _s, const VIT _e, const ST _pdc, const ST _dc, const ST _n2c)
	{
		for (VIT i = _s; i < _e; ++i)
		{
			// clear cache
			for (VIT k = _s; k < _e; ++k)
			{
				this->tempC[k][0] = 0;
				this->tempC[k][1] = 0;
				this->tempC[k][2] = 0;
			}

			for (const auto pn : vs[i][0])
			{
				for (const auto ppn : vs[pn][0])
				{
					++this->tempC[ppn][0];
					++this->tempC[ppn][1];
					++this->tempC[ppn][2];
				}
				for (const auto pnn : vs[pn][1])
				{
					++this->tempC[pnn][0];
					++this->tempC[pnn][1];
				}
			}
			for (const auto nn : vs[i][1])
			{
				for (const auto npn : vs[nn][0])
				{
					++this->tempC[npn][1];
					++this->tempC[npn][2];
				}
				for (const auto nnn : vs[nn][1])
				{
					++this->tempC[nnn][1];
				}
			}

			ST cnt = 0;
			for (VIT k = _s; k < _e; ++k)
			{
				if (k != i)
				{
					if (this->tempC[k][1] >= _dc && this->tempC[k][0] >= _pdc && this->tempC[k][2] >= _pdc)
					{
						++cnt;
					}
				}
			}

			if (cnt < _n2c)
			{
				this->tempQ.push(i);
			}
		}
	};

	init(0, vs1s, this->pdcL, this->q, this->p - 1);
	init(vs1s, va, this->pdcR, this->p, this->q - 1);

	// remove uneeded vertices, not peeling
	while (this->tempQ.size())
	{
		VIT vi = this->tempQ.front();
		this->tempQ.pop();

		for (const auto pn : vs[vi][0])
		{
			vs[pn][0].erase(vi);
		}
		for (const auto nn : vs[vi][1])
		{
			vs[nn][1].erase(vi);
		}

		vs[vi][0].clear();
		vs[vi][1].clear();
	}

	_bg.compact();
}

void mbes_opt::solve(const bigraph &_bg)
{
#ifdef mem_on
	this->ir->print_mem("graph");
#endif

	this->vs = &_bg.data();
	this->iters = 0;

	// swap pq
	if (_bg.is_swapped())
	{
		ST temp = this->p;
		this->p = this->q;
		this->q = temp;
		temp = this->pdcL;
		this->pdcL = this->pdcR;
		this->pdcR = temp;
	}
	
	this->l.clear();
	this->r.clear();
	this->cl.clear();
	this->cr.clear();
	_bg.upper(this->tempL);
	for (const auto w : this->tempL)
	{
		// the number of negative edges connected to r
		this->cl.emplace_back(vs_pair(w, 0));
	}
	_bg.lower(this->tempR);
	for (const auto w : this->tempR)
	{
		// the number of negative edges connected to l
		this->cr.emplace_back(vs_pair(w, 0));
	}
	this->l.reserve(this->cl.size());
	this->r.reserve(this->cr.size());

	pos_info pi;
	pi[0].first = 0;
	pi[0].second = (ST)this->cl.size() - 1;
	pi[1].first = 0;
	pi[1].second = (ST)this->cr.size() - 1;
	this->solve_imp(pi, 0);

	// restore pq
	if (_bg.is_swapped())
	{
		ST temp = this->p;
		this->p = this->q;
		this->q = temp;
		temp = this->pdcL;
		this->pdcL = this->pdcR;
		this->pdcR = temp;
	}

	cout << "number of iterations:=>{" << this->iters << "}" << endl;
}

ULL mbes_opt::solve_imp(const pos_info &_pInfo, const ST _status)
{
#ifdef mem_on
	this->ir->print_mem("running");
#endif
	
	ULL f = 0;

	ST ls = (ST)this->l.size(), rs = (ST)this->r.size(), np = this->p - ls, nq = this->q - rs;
	ST cls = _pInfo[0].first, cle = _pInfo[0].second, crs = _pInfo[1].first, cre = _pInfo[1].second;
	ST ccls = cle - cls + 1, ccrs = cre - crs + 1;
	pos_info npInfo;

	if (this->pf)
	{
		if (_status == 1)
		{
			if (!this->pruning_lr(this->l, this->cr, crs, cre, rs, this->pdcL))
			{
				return 0;
			}
			if (this->tempUS.size())
			{
				// move vertices from cr to r
				this->move_clcr_lr(this->cr, crs, cre, this->r);
				this->update_lr(this->tempUS, this->l);
				this->update_clcr(this->tempUS, this->cl, cls, cle);
			}

			rs = (ST)this->r.size();
			nq = this->q - rs;
			if (!this->pruning_lr(this->r, this->cl, cls, cle, ls, this->pdcR, false))
			{
				this->collect(this->cl, cls, cle);
				this->update_order(this->cl, cls, cle, cle);
				this->collect(this->cr, crs, cre);
				this->update_order(this->cr, crs, cre, cre);
				goto restore;
			}
		}
		else if (_status == 2)
		{
			if (!this->pruning_lr(this->r, this->cl, cls, cle, ls, this->pdcR))
			{
				return 0;
			}
			if (this->tempUS.size())
			{
				// move vertices from cl to l
				this->move_clcr_lr(this->cl, cls, cle, this->l);
				this->update_lr(this->tempUS, this->r);
				this->update_clcr(this->tempUS, this->cr, crs, cre);
			}

			ls = (ST)this->l.size();
			np = this->p - ls;
			if (!this->pruning_lr(this->l, this->cr, crs, cre, rs, this->pdcL, false))
			{
				this->collect(this->cl, cls, cle);
				this->update_order(this->cl, cls, cle, cle);
				this->collect(this->cr, crs, cre);
				this->update_order(this->cr, crs, cre, cre);
				goto restore;
			}
		}
		
		this->pruning_clcr(this->cl, cls, cle, this->cr, crs, cre, rs, this->pdcL);
		this->pruning_clcr(this->cr, crs, cre, this->cl, cls, cle, ls, this->pdcR);

		this->collect(this->cl, cls, cle);
		this->update_order(this->cl, cls, cle, cle);
		this->collect(this->cr, crs, cre);
		this->update_order(this->cr, crs, cre, cre);

		ccls = cle - cls + 1, ccrs = cre - crs + 1;
		if (ls + ccls < this->p || rs + ccrs < this->q)
		{
			goto restore;
		}
	}

	++this->iters;

	if (ccls || ccrs)
	{
		// special case
		if (this->sp && this->is_special(this->l, this->r, this->cl, cls, cle, this->cr, crs, cre))
		{
			this->sl.clear();
			this->sr.clear();
			this->scr.clear();
			vertex_index_list sxr;
			ST tp, tq;
			if (this->sf)
			{
				if (cre - crs + 1 <= cle - cls + 1)
				{
					for (const auto &ws : this->l)
					{
						this->sl.emplace_back(ws.first);
					}
					for (ST k = cls; k <= cle; ++k)
					{
						this->sl.emplace_back(cl[k].first);
					}
					for (const auto &ws : this->r)
					{
						this->sr.emplace_back(ws.first);
					}
					for (ST k = crs; k <= cre; ++k)
					{
						this->scr.emplace_back(cr[k].first);
					}

					tp = this->p;
					tq = this->q;
				}
				else
				{
					for (const auto &ws : this->r)
					{
						this->sl.emplace_back(ws.first);
					}
					for (ST k = crs; k <= cre; ++k)
					{
						this->sl.emplace_back(cr[k].first);
					}
					for (const auto &ws : this->l)
					{
						this->sr.emplace_back(ws.first);
					}
					for (ST k = cls; k <= cle; ++k)
					{
						this->scr.emplace_back(cl[k].first);
					}

					tp = this->q;
					tq = this->p;
				}
			}
			else
			{
				for (const auto &ws : this->l)
				{
					this->sl.emplace_back(ws.first);
				}
				for (ST k = cls; k <= cle; ++k)
				{
					this->sl.emplace_back(cl[k].first);
				}
				for (const auto &ws : this->r)
				{
					this->sr.emplace_back(ws.first);
				}
				for (ST k = crs; k <= cre; ++k)
				{
					this->scr.emplace_back(cr[k].first);
				}

				tp = this->p;
				tq = this->q;
			}

			stable_sort(this->sl.begin(), this->sl.end());

			// batch add
			ST sls = (ST)this->sl.size(), srs = (ST)this->scr.size();
			if (sls && srs)
			{
				this->intersection_special((*this->vs)[this->sl[0]], this->scr, 0, srs - 1, this->tempR);
				for (ST i = 1; i < sls; ++i)
				{
					this->intersection_special((*this->vs)[this->sl[i]], this->tempR, 0, (ST)this->tempR.size() - 1, this->tempL);
					this->tempR.swap(this->tempL);
				}

				for (const auto w : this->tempR)
				{
					this->sr.emplace_back(w);
				}
				this->difference_special(this->scr, 0, srs - 1, this->tempR, this->tempL);
				this->scr.swap(this->tempL);
			}

			ULL oa = this->idx->attempts();
			this->special(this->sl, this->sr, this->scr, sxr, tp, tq);
			if (this->idx->attempts() > oa)
			{
				f = 1;
				goto restore;
			}
			goto restore;
		}

		while (cls <= cle || crs <= cre)
		{
			if (ls + cle - cls + 1 < this->p || rs + cre - crs + 1 < this->q)
			{
				break;
			}

			VIT vil = VIT_INVALID, vir = VIT_INVALID;
			ST ncle, ncre, ts;

			if (this->sf)
			{
				if (ccls >= ccrs)
				{
					if (crs <= cre)
					{
						vir = this->cr[crs].first;
					}
					else
					{
						vil = this->cl[cls].first;
					}
				}
				else
				{
					if (cls <= cle)
					{
						vil = this->cl[cls].first;
					}
					else
					{
						vir = this->cr[crs].first;
					}
				}
			}
			else
			{
				if (crs <= cre)
				{
					vir = this->cr[crs].first;
				}
				else
				{
					vil = this->cl[cls].first;
				}
			}

			if (vil < vir)
			{
				this->l.emplace(this->cl[cls++]);
				this->intersection_clcr((*this->vs)[vil], this->cr, crs, cre);
				
				ts = 0;
				for (ST k = crs; k <= cre; ++k)
				{
					if (!masked(this->cr[k].first))
					{
						++ts;
					}
				}

				if (ts >= nq)
				{
					this->intersection_lr((*this->vs)[vil], this->r);
					this->collect(this->cr, crs, cre);
					this->update_order(this->cr, crs, cre, ncre);

					this->batch(this->cr, crs, ncre, this->cl, cls, cle);
					if (this->tempUS.size())
					{
						// move vertices from cl to l
						this->move_clcr_lr(this->cl, cls, cle, this->l);
					}
					this->collect(this->cl, cls, cle);
					this->update_order(this->cl, cls, cle, ncle);

					npInfo[0].first = cls;
					npInfo[0].second = ncle;
					npInfo[1].first = crs;
					npInfo[1].second = ncre;

					f += this->solve_imp(npInfo, 1);

					// restore r & cr
					this->intersection_lr_r((*this->vs)[vil], this->r);
					this->intersection_clcr_r((*this->vs)[vil], this->cr, crs, ncre);
					this->restore_order(this->cr, crs, ncre, ncre + 1, cre);
					// restore cl
					this->collect_ext(this->cl, ncle + 1, cle, this->tempL);
					for (const auto w : this->tempL)
					{
						this->l.erase(w);
					}
					this->restore_order(this->cl, cls, ncle, ncle + 1, cle);
				}
				else
				{
					this->collect(this->cr, crs, cre);
					this->update_order(this->cr, crs, cre, ncre);
					this->intersection_clcr_r((*this->vs)[vil], this->cr, crs, ncre);
					this->restore_order(this->cr, crs, ncre, ncre + 1, cre);
				}
				// restore l
				this->l.erase(vil);
			}
			else
			{
				this->r.emplace(this->cr[crs++]);
				this->intersection_clcr((*this->vs)[vir], this->cl, cls, cle);

				ts = 0;
				for (ST k = cls; k <= cle; ++k)
				{
					if (!masked(this->cl[k].first))
					{
						++ts;
					}
				}

				if (ts >= np)
				{
					this->intersection_lr((*this->vs)[vir], this->l);
					this->collect(this->cl, cls, cle);
					this->update_order(this->cl, cls, cle, ncle);

					this->batch(this->cl, cls, ncle, this->cr, crs, cre);
					if (this->tempUS.size())
					{
						// move vertices from cr to r
						this->move_clcr_lr(this->cr, crs, cre, this->r);
					}
					this->collect(this->cr, crs, cre);
					this->update_order(this->cr, crs, cre, ncre);

					npInfo[0].first = cls;
					npInfo[0].second = ncle;
					npInfo[1].first = crs;
					npInfo[1].second = ncre;

					f += this->solve_imp(npInfo, 2);

					// restore l & cl
					this->intersection_lr_r((*this->vs)[vir], this->l);
					this->intersection_clcr_r((*this->vs)[vir], this->cl, cls, ncle);
					this->restore_order(this->cl, cls, ncle, ncle + 1, cle);
					// restore cr
					this->collect_ext(this->cr, ncre + 1, cre, this->tempR);
					for (const auto w : this->tempR)
					{
						this->r.erase(w);
					}
					this->restore_order(this->cr, crs, ncre, ncre + 1, cre);
				}
				else
				{
					this->collect(this->cl, cls, cle);
					this->update_order(this->cl, cls, cle, ncle);
					this->intersection_clcr_r((*this->vs)[vir], this->cl, cls, ncle);
					this->restore_order(this->cl, cls, ncle, ncle + 1, cle);
				}
				// restore r
				this->r.erase(vir);
			}
		}
	}

	if (!f)
	{
		if (np <= 0 && nq <= 0)
		{
			// verify whether current biclique is a theta-biclique
			for (const auto &us : this->l)
			{
				if (rs - us.second < (ST)ceil(this->theta * rs))
				{
					f = 0;
					goto restore;
				}
			}
			for (const auto &vs : this->r)
			{
				if (ls - vs.second < (ST)ceil(this->theta * ls))
				{
					f = 0;
					goto restore;
				}
			}

			f = 1;
			this->idx->insert(this->l, this->r);
		}
		else
		{
			f = 0;
			goto restore;
		}
	}

restore:
	if (this->pf)
	{
		// restore l and r
		this->collect_ext(this->cl, cle + 1, _pInfo[0].second, this->tempL);
		this->collect_ext(this->cr, cre + 1, _pInfo[1].second, this->tempR);
		for (const auto w : this->tempL)
		{
			this->l.erase(w);
		}
		for (const auto w : this->tempR)
		{
			this->r.erase(w);
		}
		this->update_lr_r(this->tempL, this->r);
		this->update_lr_r(this->tempR, this->l);

		// restore cl and cr
		this->restore_order(this->cl, _pInfo[0].first, cle, cle + 1, _pInfo[0].second);
		this->restore_order(this->cr, _pInfo[1].first, cre, cre + 1, _pInfo[1].second);
		this->update_clcr_r(this->tempL, this->cr, _pInfo[1].first, _pInfo[1].second);
		this->update_clcr_r(this->tempR, this->cl, _pInfo[0].first, _pInfo[0].second);
	}

#ifdef mem_on
	this->ir->print_mem("running");
#endif

	return f > 0;
}

bool mbes_opt::intersection_count(const vertex_index_list &_vil, const vs_list &_vsl, ST _c2, ST _e2, const ST _target)
{
	ST e1 = (ST)_vil.size(), n2 = _e2 - _c2 + 1, pd = 0;
	if (e1 + n2 >= n2 * log2(e1))
	{
		--e1;
		for (ST k = _c2; k <= _e2; ++k)
		{
			auto &ws = _vsl[k];
			if (masked(ws.first))
			{
				continue;
			}
			if (binary_search(_vil, 0, e1, ws.first))
			{
				++pd;
				if (pd > _target)
				{
					break;
				}
			}
		}

		if (pd >= _target)
		{
			return true;
		}
		return false;
	}

	--e1;
	ST c1 = 0;
	VIT v1, v2;

	while (c1 <= e1 && _c2 <= _e2)
	{
		v1 = _vil[c1];
		v2 = _vsl[_c2].first;
		if (masked(v2))
		{
			++_c2;
		}
		else if (v1 == v2)
		{
			++pd;
			if (pd > _target)
			{
				break;
			}
			++c1;
			++_c2;
		}
		else if (v1 > v2)
		{
			++_c2;
		}
		else
		{
			++c1;
		}
	}

	if (pd >= _target)
	{
		return true;
	}
	return false;
}

ST mbes_opt::intersection_critical(const vertex_index_list &_vil, const vs_list &_vsl, ST _c2, ST _e2, const ST _target, bool _f)
{
	if (_f)
	{
		this->tempL.clear();
	}

	ST e1 = (ST)_vil.size(), n2 = _e2 - _c2 + 1, pd = 0;
	if (e1 + n2 >= n2 * log2(e1))
	{
		--e1;
		for (ST k = _c2; k <= _e2; ++k)
		{
			auto &ws = _vsl[k];
			if (masked(ws.first))
			{
				continue;
			}
			if (binary_search(_vil, 0, e1, ws.first))
			{
				++pd;
				if (_f)
				{
					this->tempL.emplace_back(ws.first);
				}
				if (pd > _target)
				{
					break;
				}
			}
		}

		return pd;
	}

	--e1;
	ST c1 = 0;
	VIT v1, v2;

	while (c1 <= e1 && _c2 <= _e2)
	{
		v1 = _vil[c1];
		v2 = _vsl[_c2].first;
		if (masked(v2))
		{
			++_c2;
		}
		else if (v1 == v2)
		{
			++pd;
			if (_f)
			{
				this->tempL.emplace_back(v1);
			}
			if (pd > _target)
			{
				break;
			}
			++c1;
			++_c2;
		}
		else if (v1 > v2)
		{
			++_c2;
		}
		else
		{
			++c1;
		}
	}

	return pd;
}

void mbes_opt::intersection_ext(const pn_neighbors &_pnn, const vs_list &_vsl, ST _c3, ST _e3, bool _inc)
{
	this->tempVS1.clear();
	const auto &pn = _pnn[0];
	const auto &nn = _pnn[1];
	ST delta = 1;
	if (!_inc)
	{
		delta = -1;
	}

	ST e1 = (ST)pn.size(), e2 = (ST)nn.size();
	ST n1 = (ST)(e1 + e2), n2 = _e3 - _c3 + 1;
	--e1;
	--e2;
	if (n1 + n2 >= n2 * log2(n1))
	{
		for (ST k = _c3; k <= _e3; ++k)
		{
			auto &ws = _vsl[k];
			if (masked(ws.first))
			{
				continue;
			}
			if (binary_search(pn, 0, e1, ws.first))
			{
				this->tempVS1.emplace_back(ws);
			}
			else if (binary_search(nn, 0, e2, ws.first))
			{
				this->tempVS1.emplace_back(vs_pair(ws.first, ws.second + delta));
			}
		}
		return;
	}

	ST c1 = 0, c2 = 0;
	VIT v1, v2, v3;
	while (c1 <= e1 && c2 <= e2)
	{
		v1 = pn[c1];
		v2 = nn[c2];
		if (v1 > v2)
		{
			while (_c3 <= _e3)
			{
				v3 = _vsl[_c3].first;
				if (masked(v3))
				{
					++_c3;
				}
				if (v3 < v2)
				{
					++_c3;
				}
				else if (v3 == v2)
				{
					this->tempVS1.emplace_back(vs_pair(v3, _vsl[_c3++].second + delta));
					break;
				}
				else
				{
					break;
				}
			}
			if (_c3 > _e3)
			{
				return;
			}
			++c2;
		}
		else
		{
			while (_c3 <= _e3)
			{
				v3 = _vsl[_c3].first;
				if (masked(v3))
				{
					++_c3;
				}
				if (v3 < v1)
				{
					++_c3;
				}
				else if (v3 == v1)
				{
					this->tempVS1.emplace_back(_vsl[_c3++]);
				}
				else
				{
					break;
				}
			}
			if (_c3 > _e3)
			{
				return;
			}
			++c1;
		}
	}

	while (c1 <= e1)
	{
		v1 = pn[c1];
		while (_c3 <= _e3)
		{
			v3 = _vsl[_c3].first;
			if (masked(v3))
			{
				++_c3;
			}
			if (v3 < v1)
			{
				++_c3;
			}
			else if (v3 == v1)
			{
				this->tempVS1.emplace_back(_vsl[_c3++]);
			}
			else
			{
				break;
			}
		}
		if (_c3 > _e3)
		{
			return;
		}
		++c1;
	}
	while (c2 <= e2)
	{
		v2 = nn[c2];
		while (_c3 <= _e3)
		{
			v3 = _vsl[_c3].first;
			if (masked(v3))
			{
				++_c3;
			}
			if (v3 < v2)
			{
				++_c3;
			}
			else if (v3 == v2)
			{
				this->tempVS1.emplace_back(vs_pair(v3, _vsl[_c3++].second+ delta));
				break;
			}
			else
			{
				break;
			}
		}
		if (_c3 > _e3)
		{
			return;
		}
		++c2;
	}
}

void mbes_opt::collect(const vs_list &_vsl, const ST _s, const ST _e)
{
	this->tempVS1.clear();
	for (ST k = _s; k <= _e; ++k)
	{
		auto &vs = _vsl[k];
		if (masked(vs.first))
		{
			this->tempVS1.emplace_back(vs_pair(unmask(vs.first), vs.second));
		}
	}
}

void mbes_opt::collect_ext(vs_list &_vsl, const ST _s, const ST _e, vertex_index_list &_con)
{
	_con.clear();
	for (ST k = _s; k <= _e; ++k)
	{
		if (masked_ext(_vsl[k].first))
		{
			_vsl[k].first = unmask_ext(_vsl[k].first);
			_con.emplace_back(_vsl[k].first);
		}
	}
}

void mbes_opt::update_order(vs_list &_vsl, ST _s, const ST _e, ST &_ne)
{
	for (ST i = _s; i <= _e; ++i)
	{
		if (!masked(_vsl[i].first))
		{
			if (_s != i)
			{
				_vsl[_s] = _vsl[i];
			}
			++_s;
		}
	}

	copy(this->tempVS1.begin(), this->tempVS1.end(), _vsl.begin() + _s);
	_ne = _s - 1;
}

void mbes_opt::restore_order(vs_list &_vsl, ST _c1, ST _e1, ST _c2, ST _e2)
{
	if (_c2 == _c1 || _c2 > _e2)
	{
		return;
	}

	this->tempVS1.clear();
	this->tempVS1.reserve((LL)_e1 - _c1 + _e2 - _c2 + 2);

	VIT v1, v2;
	ST cb1 = _c1;
	while (_c1 <= _e1 && _c2 <= _e2)
	{
		v1 = _vsl[_c1].first;
		v2 = _vsl[_c2].first;
		if (v1 > v2)
		{
			this->tempVS1.emplace_back(_vsl[_c2++]);
		}
		else
		{
			this->tempVS1.emplace_back(_vsl[_c1++]);
		}
	}
	while (_c1 <= _e1)
	{
		this->tempVS1.emplace_back(_vsl[_c1++]);
	}
	while (_c2 <= _e2)
	{
		this->tempVS1.emplace_back(_vsl[_c2++]);
	}
	copy(this->tempVS1.begin(), this->tempVS1.end(), _vsl.begin() + cb1);
}

void mbes_opt::batch(const vs_list &_cr, const ST _s1, const ST _e1, const vs_list &_cl, const ST _s2, const ST _e2)
{
	this->tempUS.clear();
	for (ST k = _s2; k <= _e2; ++k)
	{
		VIT v = _cl[k].first;
		const auto &vpn = (*this->vs)[v][0];
		ST te = (ST)vpn.size() - 1;
		bool cf = true;

		if (_cl[k].second)
		{
			continue;
		}

		for (ST i = _s1; i <= _e1; ++i)
		{
			if (!binary_search(vpn, 0, te, _cr[i].first))
			{
				cf = false;
				break;
			}
		}
		
		if (cf)
		{
			this->tempUS.emplace(v);
		}
	}
}

void mbes_opt::intersection_lr(const pn_neighbors &_pnn, unordered_map<VIT, ST> &_vsl)
{
	ST n1 = (ST)_pnn[1].size();

	if (n1 >= (ST)_vsl.size() * log2(n1))
	{
		--n1;
		for (auto &ns : _vsl)
		{
			if (binary_search(_pnn[1], 0, n1, ns.first))
			{
				++ns.second;
			}
		}
	}
	else
	{
		for (const auto nn : _pnn[1])
		{
			auto iter = _vsl.find(nn);
			if (iter != _vsl.end())
			{
				++iter->second;
			}
		}
	}
}

void mbes_opt::intersection_clcr(const pn_neighbors &_pnn, vs_list &_vsl, ST _c3, ST _e3)
{
	const auto &pn = _pnn[0];
	const auto &nn = _pnn[1];

	ST e1 = (ST)pn.size(), e2 = (ST)nn.size();
	ST n1 = (ST)(e1 + e2), n2 = _e3 - _c3 + 1;
	--e1;
	--e2;
	if (n1 + n2 >= n2 * log2(n1))
	{
		for (ST k = _c3; k <= _e3; ++k)
		{
			auto &ws = _vsl[k];
			if (binary_search(nn, 0, e2, ws.first))
			{
				++ws.second;
			}
			else if (!binary_search(pn, 0, e1, ws.first))
			{
				mask(ws.first);
			}
		}
		return;
	}

	ST c1 = 0, c2 = 0;
	VIT v1, v2, v3;
	while (c1 <= e1 && c2 <= e2)
	{
		v1 = pn[c1];
		v2 = nn[c2];
		if (v1 > v2)
		{
			while (_c3 <= _e3)
			{
				v3 = _vsl[_c3].first;
				if (v3 < v2)
				{
					mask(_vsl[_c3++].first);
				}
				else if (v3 == v2)
				{
					++_vsl[_c3++].second;
					break;
				}
				else
				{
					break;
				}
			}
			if (_c3 > _e3)
			{
				return;
			}
			++c2;
		}
		else
		{
			while (_c3 <= _e3)
			{
				v3 = _vsl[_c3].first;
				if (v3 < v1)
				{
					mask(_vsl[_c3++].first);
				}
				else if (v3 == v1)
				{
					++_c3;
				}
				else
				{
					break;
				}
			}
			if (_c3 > _e3)
			{
				return;
			}
			++c1;
		}
	}

	while (c1 <= e1)
	{
		v1 = pn[c1];
		while (_c3 <= _e3)
		{
			v3 = _vsl[_c3].first;
			if (v3 < v1)
			{
				mask(_vsl[_c3++].first);
			}
			else if (v3 == v1)
			{
				++_c3;
			}
			else
			{
				break;
			}
		}
		if (_c3 > _e3)
		{
			return;
		}
		++c1;
	}
	while (c2 <= e2)
	{
		v2 = nn[c2];
		while (_c3 <= _e3)
		{
			v3 = _vsl[_c3].first;
			if (v3 < v2)
			{
				mask(_vsl[_c3++].first);
			}
			else if (v3 == v2)
			{
				++_vsl[_c3++].second;
				break;
			}
			else
			{
				break;
			}
		}
		if (_c3 > _e3)
		{
			return;
		}
		++c2;
	}

	while (_c3 <= _e3)
	{
		mask(_vsl[_c3++].first);
	}
}

void mbes_opt::intersection_lr_r(const pn_neighbors &_pnn, unordered_map<VIT, ST> &_vsl)
{
	ST n1 = (ST)_pnn[1].size();

	if (n1 >= (ST)_vsl.size() * log2(n1))
	{
		--n1;
		for (auto &ns : _vsl)
		{
			if (binary_search(_pnn[1], 0, n1, ns.first))
			{
				--ns.second;
			}
		}
	}
	else
	{
		for (const auto nn : _pnn[1])
		{
			auto iter = _vsl.find(nn);
			if (iter != _vsl.end())
			{
				--iter->second;
			}
		}
	}
}

void mbes_opt::intersection_clcr_r(const pn_neighbors &_pnn, vs_list &_vsl, ST _c3, ST _e3)
{
	const auto &nn = _pnn[1];

	ST e1 = (ST)nn.size(), n2 = _e3 - _c3 + 1;
	if (e1 + n2 >= n2 * log2(e1))
	{
		--e1;
		for (ST k = _c3; k <= _e3; ++k)
		{
			auto &ws = _vsl[k];
			if (binary_search(nn, 0, e1, ws.first))
			{
				--ws.second;
			}
		}
		return;
	}

	--e1;
	ST c1 = 0;
	VIT v1, v3;
	while (c1 <= e1 && _c3 <= _e3)
	{
		v1 = nn[c1];
		v3 = _vsl[_c3].first;
		if (v1 == v3)
		{
			--_vsl[_c3++].second;
			++c1;
		}
		else if (v1 > v3)
		{
			++_c3;
		}
		else
		{
			++c1;
		}
	}
}

void mbes_opt::move_clcr_lr(vs_list &_src, const ST _s, const ST _e, unordered_map<VIT, ST> &_dst)
{
	for (ST k = _s; k <= _e; ++k)
	{
		auto &ws = _src[k];
		if (this->tempUS.find(ws.first) != this->tempUS.end())
		{
			_dst.emplace(ws);
			mask(ws.first);
			mask_ext(ws.first);
		}
	}
}

void mbes_opt::update_lr(const unordered_set<VIT> &_rs, unordered_map<VIT, ST> &_nl)
{
	for (const auto w : _rs)
	{
		this->intersection_lr((*this->vs)[w], _nl);
	}
}

void mbes_opt::update_clcr(const unordered_set<VIT> &_rs, vs_list &_cl, ST _c2, ST _e2)
{
	auto iter = _rs.begin();
	this->intersection_ext((*this->vs)[*iter++], _cl, _c2, _e2);
	this->tempVS2.swap(this->tempVS1);
	while (tempVS2.size() && iter != _rs.end())
	{
		this->intersection_ext((*this->vs)[*iter++], this->tempVS2, 0, (ST)this->tempVS2.size() - 1);
		this->tempVS2.swap(this->tempVS1);
	}

	ST e1 = (ST)this->tempVS2.size() - 1, c1 = 0;
	VIT v1, v2;
	while (c1 <= e1 && _c2 <= _e2)
	{
		v1 = this->tempVS2[c1].first;
		v2 = _cl[_c2].first;
		if (masked(v2))
		{
			++_c2;
		}
		else if (v1 == v2)
		{
			_cl[_c2++].second = this->tempVS2[c1].second;
			++c1;
		}
		else if (v1 > v2)
		{
			mask(_cl[_c2++].first);
		}
		else
		{
			++c1;
		}
	}

	while (_c2 <= _e2)
	{
		mask(_cl[_c2++].first);
	}
}

void mbes_opt::update_lr_r(const vertex_index_list &_rs, unordered_map<VIT, ST> &_l)
{
	for (const auto v : _rs)
	{
		this->intersection_lr_r((*this->vs)[v], _l);
	}
}

void mbes_opt::update_clcr_r(const vertex_index_list &_rs, vs_list &_cl, ST _c2, ST _e2)
{
	if (!_rs.size())
	{
		return;
	}

	auto iter = _rs.begin();
	this->intersection_ext((*this->vs)[*iter++], _cl, _c2, _e2, false);
	this->tempVS2.swap(this->tempVS1);
	while (tempVS2.size() && iter != _rs.end())
	{
		this->intersection_ext((*this->vs)[*iter++], this->tempVS2, 0, (ST)this->tempVS2.size() - 1, false);
		this->tempVS2.swap(this->tempVS1);
	}

	ST e1 = (ST)this->tempVS2.size() - 1, c1 = 0;
	VIT v1, v2;
	while (c1 <= e1 && _c2 <= _e2)
	{
		v1 = this->tempVS2[c1].first;
		v2 = _cl[_c2].first;
		if (v1 == v2)
		{
			_cl[_c2++].second = this->tempVS2[c1].second;
			++c1;
		}
		else if (v1 > v2)
		{
			++_c2;
		}
		else
		{
			++c1;
		}
	}
}

bool mbes_opt::pruning_lr(const unordered_map<VIT, ST> &_l, const vs_list &_cr, const ST _s, const ST _e, const ST _rs, const ST _pdc, bool _f)
{
	if (_f)
	{
		this->tempUS.clear();
	}

	for (const auto &ws : _l)
	{
		ST pdc = max((ST)(ws.second / this->omt), _pdc + ws.second) - _rs;

		if (pdc > 0)
		{
			ST pd = this->intersection_critical((*this->vs)[ws.first][0], _cr, _s, _e, pdc, _f);
			if (pd < pdc)
			{
				return false;
			}
			else if (_f && pd == pdc)
			{
				for (const auto vi : this->tempL)
				{
					this->tempUS.emplace(vi);
				}
			}
		}
	}

	return true;
}

void mbes_opt::pruning_clcr(vs_list &_cl, const ST _s1, const ST _e1, const vs_list &_cr, const ST _s2, const ST _e2, const ST _rs, const ST _pdc)
{
	for (ST k = _s1; k <= _e1; ++k)
	{
		auto &ws = _cl[k];
		if (masked(ws.first))
		{
			continue;
		}
		ST pdc = max((ST)(ws.second / this->omt), _pdc + ws.second) - _rs;

		if (pdc > 0 && !this->intersection_count((*this->vs)[ws.first][0], _cr, _s2, _e2, pdc))
		{
			mask(ws.first);
		}
	}
}

bool mbes_opt::is_nd(const ST _w, const vs_list &_vsl, ST _c2, ST _e2, const ST _ndc)
{
	ST nd = 0;
	const auto &v2n = (*this->vs)[_w][1];

	ST e1 = (ST)v2n.size() - 1, c1 = 0;
	VIT v1, v2;
	while (c1 <= e1 && _c2 <= _e2)
	{
		v1 = v2n[c1];
		v2 = _vsl[_c2].first;
		if (v1 == v2)
		{
			++nd;
			if (nd > _ndc)
			{
				return false;
			}
			++c1;
			++_c2;
		}
		else if (v1 > v2)
		{
			++_c2;
		}
		else
		{
			++c1;
		}
	}

	return true;
};

bool mbes_opt::is_special(const unordered_map<VIT, ST> &_l, const unordered_map<VIT, ST> &_r, const vs_list &_cl, const ST _s1, const ST _e1, const vs_list &_cr, const ST _s2, const ST _e2)
{
	ST ndcl = max((ST)_r.size(), this->q);
	ndcl -= (ST)ceil(ndcl * this->theta);
	ST ndcr = max((ST)_l.size(), this->p);
	ndcr	-= (ST)ceil(ndcr * this->theta);

	for (const auto &us : _l)
	{
		if (us.second > ndcl || !this->is_nd(us.first, _cr, _s2, _e2, ndcl - us.second))
		{
			return false;
		}
	}
	for (ST k = _s1; k <= _e1; ++k)
	{
		auto &us = _cl[k];
		if (us.second > ndcl || !this->is_nd(us.first, _cr, _s2, _e2, ndcl - us.second))
		{
			return false;
		}
	}

	for (const auto &vs : _r)
	{
		if (vs.second > ndcr || !this->is_nd(vs.first, _cl, _s1, _e1, ndcr - vs.second))
		{
			return false;
		}
	}
	for (ST k = _s2; k <= _e2; ++k)
	{
		auto &vs = _cr[k];
		if (vs.second > ndcr || !this->is_nd(vs.first, _cl, _s1, _e1, ndcr - vs.second))
		{
			return false;
		}
	}

	return true;
}

void mbes_opt::intersection_special(const pn_neighbors &_pnn, const vertex_index_list &_vil, ST _c3, ST _e3, vertex_index_list &_dst)
{
	_dst.clear();
	const auto &pn = _pnn[0];
	const auto &nn = _pnn[1];

	ST e1 = (ST)pn.size(), e2 = (ST)nn.size();
	ST n1 = (ST)(e1 + e2), n2 = _e3 - _c3 + 1;
	--e1;
	--e2;
	if (n1 + n2 >= n2 * log2(n1))
	{
		for (ST k = _c3; k <= _e3; ++k)
		{
			auto w = _vil[k];
			if (binary_search(pn, 0, e1, w) || binary_search(nn, 0, e2, w))
			{
				_dst.emplace_back(w);
			}
		}
		return;
	}

	ST c1 = 0, c2 = 0;
	VIT v1, v2, v3;
	while (c1 <= e1 && c2 <= e2)
	{
		v1 = pn[c1];
		v2 = nn[c2];
		if (v1 > v2)
		{
			while (_c3 <= _e3)
			{
				v3 = _vil[_c3];
				if (v3 < v2)
				{
					++_c3;
				}
				else if (v3 == v2)
				{
					_dst.emplace_back(v3);
					++_c3;
					break;
				}
				else
				{
					break;
				}
			}
			if (_c3 > _e3)
			{
				return;
			}
			++c2;
		}
		else
		{
			while (_c3 <= _e3)
			{
				v3 = _vil[_c3];
				if (v3 < v1)
				{
					++_c3;
				}
				else if (v3 == v1)
				{
					_dst.emplace_back(v3);
					++_c3;
					break;
				}
				else
				{
					break;
				}
			}
			if (_c3 > _e3)
			{
				return;
			}
			++c1;
		}
	}

	while (c1 <= e1)
	{
		v1 = pn[c1];
		while (_c3 <= _e3)
		{
			v3 = _vil[_c3];
			if (v3 < v1)
			{
				++_c3;
			}
			else if (v3 == v1)
			{
				_dst.emplace_back(v3);
				++_c3;
				break;
			}
			else
			{
				break;
			}
		}
		if (_c3 > _e3)
		{
			return;
		}
		++c1;
	}
	while (c2 <= e2)
	{
		v2 = nn[c2];
		while (_c3 <= _e3)
		{
			v3 = _vil[_c3];
			if (v3 < v2)
			{
				++_c3;
			}
			else if (v3 == v2)
			{
				_dst.emplace_back(v3);
				++_c3;
				break;
			}
			else
			{
				break;
			}
		}
		if (_c3 > _e3)
		{
			return;
		}
		++c2;
	}
}

void mbes_opt::difference_special(const vertex_index_list &_vil1, ST _c1, ST _e1, const vertex_index_list &_vil2, vertex_index_list &_dst)
{
	_dst.clear();

	ST e2 = (ST)_vil2.size() - 1, c2 = 0;
	VIT v1, v2;

	while (_c1 <= _e1 && c2 <= e2)
	{
		v1 = _vil1[_c1];
		v2 = _vil2[c2];
		if (v1 == v2)
		{
			++_c1;
			++c2;
		}
		else if (v1 > v2)
		{
			++c2;
		}
		else
		{
			_dst.emplace_back(v1);
			++_c1;
		}
	}

	while (_c1 <= _e1)
	{
		_dst.emplace_back(_vil1[_c1++]);
	}
}

bool mbes_opt::pruning_special(const pn_neighbors &_pnn, const vertex_index_list &_nl)
{
	const auto &pn = _pnn[0];
	const auto &nn = _pnn[1];
	if (pn.size() + nn.size() < _nl.size())
	{
		return false;
	}
	
	for (const auto v : _nl)
	{
		if (!binary_search(pn, 0, (ST)pn.size() - 1, v) && !binary_search(nn, 0, (ST)nn.size() - 1, v))
		{
			return false;
		}
	}

	return true;
}

// since it does not incur too much overhead and is similar to mbeu, the implementation is omitted.
void mbes_opt::special(const vertex_index_list &_l, vertex_index_list &_r, const vertex_index_list &_cr, vertex_index_list _xr, const ST _p, const ST _q)
{
	++this->iters;

	const ST rs = (ST)_r.size(), ls = (ST)_l.size();
	if (ls >= _p && rs >= _q)
	{
		this->idx->insert(_l, _r);
	}

	vertex_index_list nl, ncr;
	ST crs = 0, cre = (ST)_cr.size() - 1;

	while (crs <= cre)
	{
		if (rs + cre - crs + 1 < _q)
		{
			break;
		}
		
		ST ds = 1;
		VIT v = _cr[crs++];
		_r.emplace_back(v);

		// generate new l set
		this->intersection_special((*this->vs)[v], _l, 0, ls - 1, nl);
		ST nls = (ST)nl.size();

		if (nls >= _p)
		{
			// maximal test
			bool flag = false;
			for (const auto w : _xr)
			{
				flag = this->pruning_special((*this->vs)[w], nl);
				if (flag)
				{
					break;
				}
			}
			if (flag)
			{
				_r.pop_back();
				continue;
			}

			// batch add operation
			this->intersection_special((*this->vs)[nl[0]], _cr, crs, cre, this->tempR);
			for (ST i = 1; i < nls; ++i)
			{
				this->intersection_special((*this->vs)[nl[i]], this->tempR, 0, (ST)this->tempR.size() - 1, this->tempL);
				this->tempR.swap(this->tempL);
			}

			ds += (ST)this->tempR.size();
			for (const auto w : this->tempR)
			{
				_r.emplace_back(w);
			}
			this->difference_special(_cr, crs, cre, this->tempR, ncr);

			this->special(nl, _r, ncr, _xr, _p, _q);
		}
		
		for (ST i = 0; i < ds; ++i)
		{
			_r.pop_back();
		}
		_xr.emplace_back(v);
	}
}
