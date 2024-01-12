#pragma once

//#define win
#define index_on
//#define mem_on

#ifdef win
#define _CRT_SECURE_NO_WARNINGS
#endif

#ifdef win
#include <Windows.h>
#include <Psapi.h>
#else
#include <limits>
#include <unistd.h>
#include <stdio.h>
#include <string.h>
#endif
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <array>
#include <string>
#include <chrono>
#include <iostream>

#ifndef win
using std::numeric_limits;
#endif
using std::unordered_map;
using std::unordered_set;
using std::vector;
using std::array;
using std::string;
using std::chrono::system_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;
using std::cout;
using std::endl;
using std::pair;

namespace core
{
	typedef unsigned int VIT;
	typedef int ST;
	typedef long long LL;
	typedef unsigned long long ULL;

#ifdef win
	constexpr VIT VIT_INVALID= UINT_MAX;
	constexpr ST ST_INVALID = INT_MAX;
#else
	constexpr VIT VIT_INVALID = numeric_limits<VIT>::max();
	constexpr ST ST_INVALID = numeric_limits<ST>::max();
#endif

	typedef vector<VIT> vertex_index_list;
	typedef pair<VIT, ST> vs_pair;
	typedef unordered_map<VIT, ST> vs_map;
	typedef array<unordered_set<VIT>, 2> pn_neighbors;
	typedef vector<pn_neighbors> vertex_set;

	ST binary_search(const vector<ULL>&, ST, ST, ULL);
	void split(const string&, const string&, vector<string>&);
	string to_string(const vertex_index_list&);
	void core_decompose(const vector<unordered_set<VIT>>&, vector<ST>&);

	class timer
	{
	public:
		void clear();
		ULL elapse(bool = true);
		void elapse(const string&);
	private:
		system_clock::time_point lastTime, currentTime;
	};

	class info_reporter
	{
	public:
		void clear_t();
		void clear_m();
		ULL get_time();
		void print_mem(string, bool=false);
	private:
		const ULL span = 500;
#ifdef win
		const double M = 1024 * 1024;
		void get_memory_usage(string) const;
#else
		const char *data_mem = "VmRSS:";
		void print_mem(int, string) const;
#endif
		timer tt, mt;
	};
}
