#include "./bigraph.h"

using namespace core;

void bigraph::load(string _filePath)
{
	this->vs.clear();
	this->uvil.clear();
	this->lvil.clear();
	this->cache.clear();
	this->va = 0;
	this->ea = 0;
	this->pea = 0;
	this->vs1s = 0;
	this->vs2s = 0;

	ifstream ifs(_filePath);
	if (ifs.is_open())
	{
		string line;
		vector<string> temp1, temp2;

		// load bipartite graph to cache
		getline(ifs, line);
		split(line, ",", temp1);
		ST ua = atoi(temp1[0].c_str());

		while (getline(ifs, line))
		{
			split(line, ",", temp1);
			array<unordered_set<VIT>, 2> pnn;
			for (const auto &item : temp1)
			{
				split(item, "|", temp2);
				VIT vi = atoi(temp2[0].c_str());
				ST si = atoi(temp2[1].c_str());
				if (si > 0)
				{
					++this->pea;
					pnn[0].emplace(vi);
				}
				else
				{
					pnn[1].emplace(vi);
				}
			}
			this->cache.emplace_back(pnn);
			++this->va;
			this->ea += (ST)temp1.size();
		}
		this->ea >>= 1;
		this->pea >>= 1;
		ifs.close();

		this->vs1s = ua;
		this->vs2s = this->va - ua;
	}
}

void bigraph::print_info() const
{
	cout << "v1:" << this->vs1s << endl;
	cout << "v2:" << this->vs2s << endl;
	cout << "e:" << this->ea << endl;
	cout << "pe:" << this->pea << endl;
	cout << "ne:" << (this->ea - this->pea) << endl;
	cout << "d:" << 2.0 * this->ea / this->va << endl;
}

void bigraph::reorder(sort_method _sm)
{
	ST va = this->va, vs1s = this->vs1s;
	auto &ovs = this->cache;
	this->vs.resize(va);

	// sort
	vs_list vsl;
	this->sort_p(_sm, vsl);

	// reorder
	unordered_map<VIT, VIT> om;
	for (ST i = 0; i < va; ++i)
	{
		om.emplace(vsl[i].first, (VIT)i);
	}

	for (ST i = 0; i < va; ++i)
	{
		auto &nv = this->vs[om[i]];
		for (const auto v : ovs[i][0])
		{
			nv[0].emplace_back(om[v]);
		}
		for (const auto &v : ovs[i][1])
		{
			nv[1].emplace_back(om[v]);
		}

		stable_sort(nv[0].begin(), nv[0].end());
		stable_sort(nv[1].begin(), nv[1].end());
	}

	// release unneeded memory
	ovs.clear();
	ovs.shrink_to_fit();

	// set upper and lower
	if (vs1s >= this->vs2s)
	{
		for (ST i = 0; i < vs1s; ++i)
		{
			this->uvil.emplace_back(i);
		}
		for (ST i = vs1s; i < va; ++i)
		{
			this->lvil.emplace_back(i);
		}
	}
	else
	{
		for (ST i = 0; i < vs1s; ++i)
		{
			this->lvil.emplace_back(i);
		}
		for (ST i = vs1s; i < va; ++i)
		{
			this->uvil.emplace_back(i);
		}
	}
}

void bigraph::compact()
{
	auto &vs = this->cache;
	ST cur = 0, nvs1s = 0, n = this->va, vs1s = this->vs1s;
	unordered_map<VIT, VIT> om;
	for (ST i = 0; i < n; ++i)
	{
		if (i == vs1s)
		{
			nvs1s = cur;
		}

		// must ensure that before calling this function, the neighbors of unneeded vertices are removed 
		if (vs[i][0].size() + vs[i][1].size())
		{
			if (i != cur)
			{
				vs[cur].swap(vs[i]);
			}
			om.emplace(i, cur);
			++cur;
		}
	}

	n = this->va = cur;
	this->vs1s = nvs1s;
	this->vs2s = cur - nvs1s;
	vs.erase(vs.begin() + cur, vs.end());

	ST ea = 0, pea = 0;
	for (ST i = 0; i < n; ++i)
	{
		auto &v = vs[i];
		array<unordered_set<VIT>, 2> pnn;
		ea += (ST)(v[0].size() + v[1].size());
		for (const auto vi : v[0])
		{
			pnn[0].emplace(om[vi]);
			++pea;
		}
		for (const auto vi : v[1])
		{
			pnn[1].emplace(om[vi]);
		}
		v.swap(pnn);
	}
	this->ea = ea >> 1;
	this->pea = pea >> 1;
}

void bigraph::write(string _path) const
{
	ofstream ofs;
	ofs.open(_path);
	if (ofs.is_open())
	{
		ofs << this->vs1s << "," << this->vs2s << endl;
		for (ST i = 0; i < this->va; ++i)
		{
			bool f = false;
			for (const auto pn : this->vs[i][0])
			{
				if (!f)
				{
					f = true;
				}
				else
				{
					ofs << ",";
				}
				ofs << pn << "|" << 1;
			}
			for (const auto nn : this->vs[i][1])
			{
				if (!f)
				{
					f = true;
				}
				else
				{
					ofs << ",";
				}
				ofs << nn << "|" << 0;
			}
			ofs << endl;
		}
	}
}

void bigraph::sort_p(sort_method _sm, vs_list &_vp_list)
{
	ST va = this->va, vs1s = this->vs1s;
	auto &vs = this->cache;
	_vp_list.resize(va);

	if (_sm == sort_method::degree_o)
	{
		for (ST i = 0; i < va; ++i)
		{
			_vp_list[i] = { i, (ST)(vs[i][0].size() + vs[i][1].size())};
		}
	}
	else if (_sm == sort_method::degree_positive_o)
	{
		for (ST i = 0; i < va; ++i)
		{
			_vp_list[i] = { i, (ST)vs[i][0].size()};
		}
	}
	else if (_sm == sort_method::degeneracy_o)
	{
		vector<unordered_set<VIT>> usg(va);
		vector<ST> cn(va);
		for (ST i = 0; i < va; ++i)
		{
			for (const auto v : vs[i][0])
			{
				usg[i].emplace(v);
			}
			for (const auto v : vs[i][1])
			{
				usg[i].emplace(v);
			}
		}

		core_decompose(usg, cn);
		for (ST i = 0; i < va; ++i)
		{
			_vp_list[i] = { i, (ST)cn[i] };
		}
	}
	else if (_sm == sort_method::degeneracy_positive_o)
	{
		vector<unordered_set<VIT>> usg(va);
		vector<ST> cn(va);
		for (ST i = 0; i < va; ++i)
		{
			for (const auto v : vs[i][0])
			{
				usg[i].emplace(v);
			}
		}

		core_decompose(usg, cn);
		for (ST i = 0; i < va; ++i)
		{
			_vp_list[i] = { i, (ST)cn[i] };
		}
	}
	else
	{
		for (ST i = 0; i < vs1s; ++i)
		{
			_vp_list[i] = { i, (ST)i };
		}
		for (ST i = vs1s; i < va; ++i)
		{
			_vp_list[i] = { i, (ST)(i - vs1s) };
		}
	}

	stable_sort(_vp_list.begin(), _vp_list.begin() + vs1s, vp_sort);
	stable_sort(_vp_list.begin() + vs1s, _vp_list.end(), vp_sort);
}
