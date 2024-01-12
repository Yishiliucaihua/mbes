#include "./utils.h"

using namespace core;

ST core::binary_search(const vector<ULL> &_arr, ST _s, ST _e, ULL _key)
{
	ST mid, ret = _e + 1;
	while (_s <= _e)
	{
		mid = (_s + _e) >> 1;
		if (_arr[mid] == _key)
		{
			return mid;
		}
		else if (_arr[mid] > _key)
		{
			_e = mid - 1;
			ret = mid;
		}
		else
		{
			_s = mid + 1;
		}
	}
	return ret;
}

void core::split(const string &_str, const string &_del, vector<string> &_container)
{
	_container.clear();
	if (_str.size())
	{
		char *strs = new char[_str.length() + 1];
		strcpy(strs, _str.c_str());

		char *d = new char[_del.length() + 1];
		strcpy(d, _del.c_str());

		char *p = strtok(strs, d);
		while (p)
		{
			_container.emplace_back(string(p));
			p = strtok(NULL, d);
		}

		delete[] strs;
		delete[] d;
	}
}

static ST bs = 1025;
static char *buffer = new char[1025];
string core::to_string(const vertex_index_list &_vil)
{
	ST n = (ST)_vil.size();
	ST max = n * 11;

	if (bs < max)
	{
		delete[] buffer;
		bs = max;
		buffer = new char[max];
	}

	char *pb = buffer;
	memset(pb, 0, max);
	ST len = snprintf(pb, max, "%u", _vil[0]);
	pb += len;
	max -= len;
	for (ST i = 1; i < n; ++i)
	{
		len = snprintf(pb, max, "_%u", _vil[i]);
		pb += len;
		max -= len;
	}

	return string(buffer);
}

void core::core_decompose(const vector<unordered_set<VIT>> &_sg, vector<ST> &_deg)
{
	ST n = (ST)_sg.size();
	vector<VIT> vert(n), pos(n);
	vector<VIT> bin;

	// init deg
	ST md = 0;
	for (ST i = 0; i < n; ++i)
	{
		_deg[i] = (ST)_sg[i].size();
		if (_deg[i] > md)
		{
			md = _deg[i];
		}
	}

	bin.resize((ULL)md + 1, 0);

	// insert to bins
	for (ST i = 0; i < n; ++i)
	{
		++bin[_deg[i]];
	}

	// set start pos of each bin
	ST start = 0, num;
	for (ST i = 0; i <= md; ++i)
	{
		num = bin[i];
		bin[i] = start;
		start += num;
	}

	// associate bin, pos and vert
	for (ST i = 0; i < n; ++i)
	{
		pos[i] = bin[_deg[i]]++;
		vert[pos[i]] = i;
	}

	// recover bin
	for (ST i = md; i >= 1; --i)
	{
		bin[i] = bin[(ULL)i - 1];
	}
	bin[0] = 0;

	// peeling
	VIT v, w;
	ST du, pu, pw;
	for (ST i = 0; i < n; ++i)
	{
		v = vert[i];
		for (const auto u : _sg[v])
		{
			if (_deg[u] > _deg[v])
			{
				du = _deg[u];
				pu = pos[u];
				pw = bin[du];
				w = vert[pw];

				if (u != w)
				{
					pos[u] = pw;
					vert[pu] = w;
					pos[w] = pu;
					vert[pw] = u;
				}

				++bin[du];
				--_deg[u];
			}
		}
	}
}

void timer::clear()
{
	this->currentTime = this->lastTime = system_clock::now();
}

ULL timer::elapse(bool _clear)
{
	this->currentTime = system_clock::now();
	ULL ret = duration_cast<milliseconds>(this->currentTime - this->lastTime).count();
	if (_clear)
	{
		this->lastTime = this->currentTime;
	}
	return ret;
}

void timer::elapse(const string &_str)
{
	ULL elapse = this->elapse();
	cout << _str << " elapse time: " << elapse << endl;
}

void info_reporter::clear_t()
{
	this->tt.clear();
}

void info_reporter::clear_m()
{
	this->mt.clear();
}

ULL info_reporter::get_time()
{
	return this->tt.elapse();
}

void info_reporter::print_mem(string _prefix, bool _force)
{
	if (_force || this->mt.elapse(false) > span)
	{
#ifdef win
		get_memory_usage(_prefix);
#else
		print_mem(getpid(), _prefix);
#endif
		this->mt.clear();
	}
}

#ifdef win
void info_reporter::get_memory_usage(string _prefix) const
{
	PROCESS_MEMORY_COUNTERS pmc;
	if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
	{
		uint64_t ms = pmc.PagefileUsage;
		cout << _prefix << " mem:=>{" << ms << "}" << endl;
	}
}
#else
void info_reporter::print_mem(int _pid, string _prefix) const
{
	FILE *stream;
	char cache[512];
	char mem_info[64];

	sprintf(mem_info, "/proc/%d/status", _pid);
	stream = fopen(mem_info, "r");
	if (stream == NULL)
	{
		return;
	}

	while (fscanf(stream, "%s", cache) != EOF)
	{
		if (strncmp(cache, data_mem, strlen(data_mem)) == 0)
		{
			if (fscanf(stream, "%s", cache) != EOF)
			{
				cout << _prefix << " mem:=>{" << string(cache) << "}" << endl;
				break;
			}
		}
	}
}
#endif
