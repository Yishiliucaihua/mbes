#include "./mbes.h"

using namespace core;

void mbes::reduction_1hop(bigraph &_bg)
{
	// clear queue
	while (this->tempQ.size())
	{
		this->tempQ.pop();
	}

	// initialize queue
	const VIT va = _bg.va, vs1s = _bg.vs1s;
	auto &vs = _bg.vs;

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

void mbes::reduction_2hop(bigraph &_bg)
{
	// clear queue
	while (this->tempQ.size())
	{
		this->tempQ.pop();
	}

	// initialize queue
	const VIT va = _bg.va, vs1s = _bg.vs1s;
	auto &vs = _bg.vs;
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

void mbes::solve(const bigraph &_bg)
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

	vs_map l, r, cl, cr;
	_bg.upper(this->tempV);
	for (const auto w : this->tempV)
	{
		// the number of negative edges connected to r
		cl.emplace(w, 0);
	}
	_bg.lower(this->tempV);
	for (const auto w : this->tempV)
	{
		// the number of negative edges connected to l
		cr.emplace(w, 0);
	}

	this->solve_imp(l, r, cl, cr);

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

ULL mbes::solve_imp(const vs_map &_l, const vs_map &_r, const vs_map &_cl, const vs_map &_cr)
{
#ifdef mem_on
	this->ir->print_mem("running");
#endif

	ULL f = 0;

	const ST ls = (ST)_l.size(), rs = (ST)_r.size();
	const ST np = this->p - ls, nq = this->q - rs;
	vs_map cl, cr;

	if (this->pf)
	{
		this->pruning_clcr(_cl, rs, _cr, this->pdcL);
		for (const auto &us : _cl)
		{
			if (this->tempUS.find(us.first) == this->tempUS.end())
			{
				cl.emplace(us);
			}
		}
		this->pruning_clcr(_cr, ls, cl, this->pdcR);
		for (const auto &vs : _cr)
		{
			if (this->tempUS.find(vs.first) == this->tempUS.end())
			{
				cr.emplace(vs);
			}
		}

		if ((ST)(ls + cl.size()) < this->p || (ST)(rs + cr.size()) < this->q)
		{
			return 0;
		}
	}
	else
	{
		cl = _cl;
		cr = _cr;
	}

	++this->iters;

	if (cl.size() || cr.size())
	{
		// special case
		if (this->sp && this->is_special(_l, _r, cl, cr))
		{
			this->tl.clear();
			this->tr.clear();
			this->tcr.clear();
			ST tp, tq;

			if (this->sf)
			{
				if (cr.size() <= cl.size())
				{
					this->tl.reserve(cl.size() + _l.size());
					this->tr.reserve(_r.size());
					this->tcr.reserve(cr.size());

					for (const auto &ws : _l)
					{
						this->tl.emplace(ws.first);
					}
					for (const auto &ws : cl)
					{
						this->tl.emplace(ws.first);
					}
					for (const auto &ws : _r)
					{
						this->tr.emplace(ws.first);
					}
					for (const auto &ws : cr)
					{
						this->tcr.emplace(ws.first);
					}

					tp = this->p;
					tq = this->q;
				}
				else
				{
					this->tl.reserve(cr.size() + _r.size());
					this->tr.reserve(_l.size());
					this->tcr.reserve(cl.size());

					for (const auto &ws : _r)
					{
						this->tl.emplace(ws.first);
					}
					for (const auto &ws : cr)
					{
						this->tl.emplace(ws.first);
					}
					for (const auto &ws : _l)
					{
						this->tr.emplace(ws.first);
					}
					for (const auto &ws : cl)
					{
						this->tcr.emplace(ws.first);
					}

					tp = this->q;
					tq = this->p;
				}
			}
			else
			{
				this->tl.reserve(cl.size() + _l.size());
				this->tr.reserve(_r.size());
				this->tcr.reserve(cr.size());

				for (const auto &ws : _l)
				{
					this->tl.emplace(ws.first);
				}
				for (const auto &ws : cl)
				{
					this->tl.emplace(ws.first);
				}
				for (const auto &ws : _r)
				{
					this->tr.emplace(ws.first);
				}
				for (const auto &ws : cr)
				{
					this->tcr.emplace(ws.first);
				}

				tp = this->p;
				tq = this->q;
			}

			// batch add
			this->tempS1 = this->tcr;
			for (const auto w1 : this->tl)
			{
				const auto &pn = (*this->vs)[w1][0];
				const auto &nn = (*this->vs)[w1][1];
				this->tempS2.clear();
				for (const auto w2 : this->tempS1)
				{
					if (pn.find(w2) != pn.end() || nn.find(w2) != nn.end())
					{
						this->tempS2.emplace(w2);
					}
				}
				this->tempS1.swap(this->tempS2);
			}
			for (const auto w : this->tempS1)
			{
				this->tr.emplace(w);
				this->tcr.erase(w);
			}

			ULL oa = this->idx->attempts();
			this->special(this->tl, this->tr, this->tcr, this->txr, tp, tq);
			if (this->idx->attempts() > oa)
			{
				return 1;
			}
			return 0;
		}

		// sort candidates according to the priority
		// we use quick sort rather than bin sort, because it is nearly linear for small vector
		list<VIT> priL, priR;
		this->tempV.clear();
		for (const auto &us : cl)
		{
			this->tempV.emplace_back(us.first);
		}
		stable_sort(this->tempV.begin(), this->tempV.end());
		for (const auto u : this->tempV)
		{
			priL.emplace_back(u);
		}

		this->tempV.clear();
		for (const auto &vs : cr)
		{
			this->tempV.emplace_back(vs.first);
		}
		stable_sort(this->tempV.begin(), this->tempV.end());
		for (const auto v : this->tempV)
		{
			priR.emplace_back(v);
		}

		const ST cls = (ST)priL.size(), crs = (ST)priR.size();
		vs_map nl, nr, ncl, ncr;

		while (priL.size() || priR.size())
		{
			if ((ST)(ls + priL.size()) < this->p || (ST)(rs + priR.size()) < this->q)
			{
				break;
			}

			VIT vil = VIT_INVALID, vir = VIT_INVALID;
			
			if (this->sf)
			{
				if (cls >= crs)
				{
					if (priR.size())
					{
						vir = priR.front();
					}
					else
					{
						vil = priL.front();
					}
				}
				else
				{
					if (priL.size())
					{
						vil = priL.front();
					}
					else
					{
						vir = priR.front();
					}
				}
			}
			else
			{
				if (priR.size())
				{
					vir = priR.front();
				}
				else
				{
					vil = priL.front();
				}
			}

			if (vil < vir)
			{
				nl = _l;
				nl.emplace(vil, cl[vil]);
				priL.pop_front();
				cl.erase(vil);

				this->intersection_clcr((*this->vs)[vil], cr, ncr);

				if ((ST)ncr.size() >= nq)
				{
					nr = _r;
					this->intersection_lr((*this->vs)[vil], nr);

					ncl = cl;
					this->tempUS.clear();
					for (const auto &v : ncl)
					{
						if (v.second)
						{
							continue;
						}

						bool cf = true;
						const auto &vpn = (*this->vs)[v.first][0];
						for (const auto &w : ncr)
						{
							if (vpn.find(w.first) == vpn.end())
							{
								cf = false;
								break;
							}
						}
						
						if (cf)
						{
							this->tempUS.emplace(v.first);
						}
					}

					for (const auto w : this->tempUS)
					{
						nl.emplace(w, ncl[w]);
						ncl.erase(w);
					}

					if (this->pf)
					{
						if (!this->pruning_lr(nl, rs, ncr, this->pdcL))
						{
							continue;
						}
						if (this->tempUS.size())
						{
							this->update_clcr(this->tempUS, ncl, this->tempVM1);
							ncl.swap(this->tempVM1);
							this->update_lr(this->tempUS, nl);
							for (const auto w : this->tempUS)
							{
								nr.emplace(w, ncr[w]);
								ncr.erase(w);
							}
						}

						if (!this->pruning_lr(nr, (ST)nl.size(), ncl, this->pdcR))
						{
							continue;
						}
					}

					f += this->solve_imp(nl, nr, ncl, ncr);
				}
			}
			else
			{
				nr = _r;
				nr.emplace(vir, cr[vir]);
				priR.pop_front();
				cr.erase(vir);

				this->intersection_clcr((*this->vs)[vir], cl, ncl);

				if ((ST)ncl.size() >= np)
				{
					nl = _l;
					this->intersection_lr((*this->vs)[vir], nl);

					ncr = cr;
					this->tempUS.clear();
					for (const auto &v : ncr)
					{
						if (v.second)
						{
							continue;
						}

						bool cf = true;
						const auto &vpn = (*this->vs)[v.first][0];
						for (const auto &w : ncl)
						{
							if (vpn.find(w.first) == vpn.end())
							{
								cf = false;
								break;
							}
						}
						
						if (cf)
						{
							this->tempUS.emplace(v.first);
						}
					}

					for (const auto w : this->tempUS)
					{
						nr.emplace(w, ncr[w]);
						ncr.erase(w);
					}

					if (this->pf)
					{
						if (!this->pruning_lr(nr, ls, ncl, this->pdcR))
						{
							continue;
						}
						if (this->tempUS.size())
						{
							this->update_clcr(this->tempUS, ncr, this->tempVM1);
							ncr.swap(this->tempVM1);
							this->update_lr(this->tempUS, nr);
							for (const auto w : this->tempUS)
							{
								nl.emplace(w, ncl[w]);
								ncl.erase(w);
							}
						}

						if (!this->pruning_lr(nl, (ST)nr.size(), ncr, this->pdcL))
						{
							continue;
						}
					}

					f += this->solve_imp(nl, nr, ncl, ncr);
				}
			}
		}
	}

	if (!f)
	{
		if (np <= 0 && nq <= 0)
		{
			// verify whether current biclique is a theta-biclique
			for (const auto &us : _l)
			{
				if (rs - us.second < (ST)ceil(this->theta * rs))
				{
					return 0;
				}
			}
			for (const auto &vs : _r)
			{
				if (ls - vs.second < (ST)ceil(this->theta * ls))
				{
					return 0;
				}
			}

			this->idx->insert(_l, _r);
		}
		else
		{
			return 0;
		}
	}

#ifdef mem_on
	this->ir->print_mem("running");
#endif

	return 1;
}

void mbes::intersection_lr(const pn_neighbors &_pnn, vs_map &_vm)
{
	const auto &nn = _pnn[1];
	if (_vm.size() <= nn.size())
	{
		for (auto &ws : _vm)
		{
			if (nn.find(ws.first) != nn.end())
			{
				++ws.second;
			}
		}
	}
	else
	{
		for (const auto w : nn)
		{
			auto wi = _vm.find(w);
			if (wi != _vm.end())
			{
				++wi->second;
			}
		}
	}
}

void mbes::intersection_clcr(const pn_neighbors &_pnn, const vs_map &_vm, vs_map &_nvm)
{
	_nvm.clear();
	const auto &pn = _pnn[0];
	const auto &nn = _pnn[1];

	if (_vm.size() <= pn.size() + nn.size())
	{
		for (const auto &ws : _vm)
		{
			if (pn.find(ws.first) != pn.end())
			{
				_nvm.emplace(ws);
			}
			else if (nn.find(ws.first) != nn.end())
			{
				_nvm.emplace(ws.first, ws.second + 1);
			}
		}
	}
	else
	{
		for (const auto w : pn)
		{
			auto wi = _vm.find(w);
			if (wi != _vm.end())
			{
				_nvm.emplace(w, wi->second);
			}
		}
		for (const auto w : nn)
		{
			auto wi = _vm.find(w);
			if (wi != _vm.end())
			{
				_nvm.emplace(w, wi->second + 1);
			}
		}
	}
}

void mbes::update_lr(const unordered_set<VIT> &_rs, vs_map &_nl)
{
	for (const auto w : _rs)
	{
		this->intersection_lr((*this->vs)[w], _nl);
	}
}

void mbes::update_clcr(const unordered_set<VIT> &_rs, const vs_map &_cl, vs_map &_ncl)
{
	auto iter = _rs.begin();
	this->intersection_clcr((*this->vs)[*iter++], _cl, _ncl);
	while (iter != _rs.end())
	{
		this->intersection_clcr((*this->vs)[*iter++], _ncl, this->tempVM2);
		_ncl.swap(this->tempVM2);
	}
}

bool mbes::pruning_lr(const vs_map &_l, const ST _rs, const vs_map &_cr, const ST _pdc)
{
	this->tempUS.clear();

	for (const auto &ws : _l)
	{
		ST pdc = max((ST)(ws.second / this->omt), _pdc + ws.second) - _rs;

		if (pdc > 0)
		{
			ST pd = 0;
			this->tempV.clear();
			const auto &v2p = (*this->vs)[ws.first][0];

			if (v2p.size() <= _cr.size())
			{
				for (const auto pn : v2p)
				{
					if (_cr.find(pn) != _cr.end())
					{
						++pd;
						this->tempV.emplace_back(pn);
						if (pd > pdc)
						{
							break;
						}
					}
				}
			}
			else
			{
				for (const auto &vs : _cr)
				{
					if (v2p.find(vs.first) != v2p.end())
					{
						++pd;
						this->tempV.emplace_back(vs.first);
						if (pd > pdc)
						{
							break;
						}
					}
				}
			}
			if (pd < pdc)
			{
				return false;
			}
			else if (pd == pdc)
			{
				for (const auto vi : this->tempV)
				{
					this->tempUS.emplace(vi);
				}
			}
		}
	}

	return true;
}

void mbes::pruning_clcr(const vs_map &_cl, const ST _rs, const vs_map &_cr, const ST _pdc)
{
	this->tempUS.clear();

	for (const auto &ws : _cl)
	{
		ST pdc = max((ST)(ws.second / this->omt), _pdc + ws.second) - _rs;

		if (pdc > 0)
		{
			ST pd = 0;
			const auto &v2p = (*this->vs)[ws.first][0];

			if (v2p.size() <= _cr.size())
			{
				for (const auto pn : v2p)
				{
					if (_cr.find(pn) != _cr.end())
					{
						++pd;
						if (pd > pdc)
						{
							break;
						}
					}
				}
			}
			else
			{
				for (const auto &vs : _cr)
				{
					if (v2p.find(vs.first) != v2p.end())
					{
						++pd;
						if (pd > pdc)
						{
							break;
						}
					}
				}
			}
			if (pd < pdc)
			{
				this->tempUS.emplace(ws.first);
			}
		}
	}
}

bool mbes::is_nd(const ST _w, const vs_map &_vm, const ST _ndc)
{
	ST nd = 0;
	const auto &v2n = (*this->vs)[_w][1];

	if (v2n.size() <= _vm.size())
	{
		for (const auto nn : v2n)
		{
			if (_vm.find(nn) != _vm.end())
			{
				++nd;
				if (nd > _ndc)
				{
					return false;
				}
			}
		}
	}
	else
	{
		for (const auto &vs : _vm)
		{
			if (v2n.find(vs.first) != v2n.end())
			{
				++nd;
				if (nd > _ndc)
				{
					return false;
				}
			}
		}
	}
	return true;
};

bool mbes::is_special(const vs_map &_l, const vs_map &_r, const vs_map &_cl, const vs_map &_cr)
{
	ST ndcl = max((ST)_r.size(), this->q);
	ndcl -= (ST)ceil(ndcl * this->theta);
	ST ndcr = max((ST)_l.size(), this->p);
	ndcr	-= (ST)ceil(ndcr * this->theta);

	for (const auto &us : _l)
	{
		if (us.second > ndcl || !this->is_nd(us.first, _cr, ndcl - us.second))
		{
			return false;
		}
	}
	for (const auto &us : _cl)
	{
		if (us.second > ndcl || !this->is_nd(us.first, _cr, ndcl - us.second))
		{
			return false;
		}
	}

	for (const auto &vs : _r)
	{
		if (vs.second > ndcr || !this->is_nd(vs.first, _cl, ndcr - vs.second))
		{
			return false;
		}
	}
	for (const auto &vs : _cr)
	{
		if (vs.second > ndcr || !this->is_nd(vs.first, _cl, ndcr - vs.second))
		{
			return false;
		}
	}

	return true;
}

void mbes::special(const unordered_set<VIT> &_l, const unordered_set<VIT> &_r, const unordered_set<VIT> &_cr, vertex_index_list _xr, const ST _p, const ST _q)
{
	++this->iters;

	const ST rs = (ST)_r.size();
	if ((ST)_l.size() >= _p && rs >= _q)
	{
		this->idx->insert(_l, _r);
	}

	unordered_set<VIT> cr = _cr;
	list<VIT> priR;
	this->tempV.clear();
	for (const auto v : cr)
	{
		this->tempV.emplace_back(v);
	}
	stable_sort(this->tempV.begin(), this->tempV.end());
	for (const auto v : this->tempV)
	{
		priR.emplace_back(v);
	}

	unordered_set<VIT> nl, nr, ncr;

	while (priR.size())
	{
		if ((ST)(rs + priR.size()) < _q)
		{
			break;
		}

		VIT v = priR.front();
		priR.pop_front();
		cr.erase(v);
		nr = _r;
		nr.emplace(v);
		ncr = cr;

		// generate new l set
		nl.clear();
		const auto &v2p = (*this->vs)[v][0];
		const auto &v2n = (*this->vs)[v][1];
		for (const auto w : _l)
		{
			if (v2p.find(w) != v2p.end() || v2n.find(w) != v2n.end())
			{
				nl.emplace(w);
			}
		}

		if ((ST)nl.size() >= _p)
		{
			// maximal test
			bool flag = false;
			for (const auto w1 : _xr)
			{
				const auto &pn = (*this->vs)[w1][0];
				const auto &nn = (*this->vs)[w1][1];
				flag = true;
				for (const auto w2 : nl)
				{
					if (pn.find(w2) == pn.end() && nn.find(w2) == nn.end())
					{
						flag = false;
						break;
					}
				}
				if (flag)
				{
					break;
				}
			}
			if (flag)
			{
				continue;
			}

			// batch add operation
			this->tempS1 = ncr;
			for (const auto w1 : nl)
			{
				const auto &pn = (*this->vs)[w1][0];
				const auto &nn = (*this->vs)[w1][1];
				this->tempS2.clear();
				for (const auto w2 : this->tempS1)
				{
					if (pn.find(w2) != pn.end() || nn.find(w2) != nn.end())
					{
						this->tempS2.emplace(w2);
					}
				}
				this->tempS1.swap(this->tempS2);
			}
			for (const auto w : this->tempS1)
			{
				nr.emplace(w);
				ncr.erase(w);
			}

			this->special(nl, nr, ncr, _xr, _p, _q);
		}

		_xr.emplace_back(v);
	}
}
