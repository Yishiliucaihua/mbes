#include "./index.h"

using namespace core;

void index_opt::init(ST _size)
{
	this->data.resize(_size);
	for (auto &item: this->data)
	{
		item.clear();
	}
	this->vil_data.clear();
	this->s = 0;
	this->iters = 0;
}

void index_opt::insert(const vs_map &_l, const vs_map &_r)
{
	this->temp.clear();
	for (const auto &v : _l)
	{
		this->temp.emplace_back(v.first);
	}
	for (const auto &v : _r)
	{
		this->temp.emplace_back(v.first);
	}

	stable_sort(this->temp.begin(), this->temp.end());
	this->add(this->temp);
}

void index_opt::insert(const unordered_set<VIT> &_l, const unordered_set<VIT> &_r)
{
	this->temp.clear();
	for (const auto v : _l)
	{
		this->temp.emplace_back(v);
	}
	for (const auto v : _r)
	{
		this->temp.emplace_back(v);
	}

	stable_sort(this->temp.begin(), this->temp.end());
	this->add(this->temp);
}

bool index_opt::insert(const vertex_index_list &_vil)
{
	const ST n = (ST)_vil.size();
	if (!n)
	{
		return false;
	}
	else if (n == 1)
	{
		if (this->data[_vil[0]].size())
		{
			return false;
		}
		else
		{
			this->data[_vil[0]].emplace_back(this->s++);
			return true;
		}
	}

	vector<const vector<ULL>*> heads;
	vector<ST> curs, sizes;
	for (const auto vi : _vil)
	{
		heads.emplace_back(&this->data[vi]);
		curs.emplace_back(0);
		sizes.emplace_back((ST)this->data[vi].size());
	}

	LL maxV = -1;
	ST maxI = -1;
	bool isDuplicate = true;
	for (ST k = 0; k < n; ++k)
	{
		if (!sizes[k])
		{
			isDuplicate = false;
			break;
		}
		if ((LL)(*heads[k])[0] > maxV)
		{
			maxV = (*heads[k])[0];
			maxI = k;
		}
	}
	while (isDuplicate)
	{
		bool f = true;
		for (ST k = 0; k < n; ++k)
		{
			if (k != maxI)
			{
				const auto &list = *heads[k];
				if ((LL)list[curs[k]] != maxV)
				{
					f = false;
					curs[k] = binary_search(list, curs[k], sizes[k] - 1, maxV);
					if (curs[k] >= sizes[k])
					{
						isDuplicate = false;
						f = true;
						break;
					}
				}
			}
		}
		if (f)
		{
			break;
		}

		maxV = -1;
		for (ST k = 0; k < n; ++k)
		{
			if ((LL)(*heads[k])[curs[k]] > maxV)
			{
				maxV = (*heads[k])[curs[k]];
				maxI = k;
			}
		}
	}
	if (isDuplicate)
	{
		return false;
	}

	for (const auto vi : _vil)
	{
		this->data[vi].emplace_back(this->s);
	}
	++this->s;
	return true;
}

void index_opt::add(const vertex_index_list &_vil)
{
	++this->iters;

#ifndef index_on
	return;
#endif

	this->vil_data.emplace_back(_vil);
}

void index_opt::remove_duplicate()
{
	for (const auto &vil : this->vil_data)
	{
		this->insert(vil);
	}
}
