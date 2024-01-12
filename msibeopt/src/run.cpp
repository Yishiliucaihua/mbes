#include "./mbes_opt.h"

using core::ST;
using core::ULL;
using core::bigraph;
using core::index_opt;
using core::info_reporter;
using core::search_interface;
using core::mbes_opt;

constexpr int parameter_num = 10;

int main(int _argc, char *_argv[])
{
	if (_argc < parameter_num)
	{
		cout << "error: miss necessary parameter" << endl;
		exit(-1);
	}
	
	info_reporter ir;
	ULL reduction_time = 0, reorder_time = 0, run_time = 0, idx_time = 0;

	ST p = atoi(_argv[2]), q = atoi(_argv[3]);
	double theta = strtod(_argv[4], nullptr);
	ST rd = atoi(_argv[5]), sf = atoi(_argv[6]), pf = atoi(_argv[7]), sp = atoi(_argv[8]);
	bigraph::sort_method sm = (bigraph::sort_method)strtoul(_argv[9], nullptr, 10);

	// load
	bigraph bg;
	bg.load(_argv[1]);
	cout << "original graph:=>{";
	bg.print_info();
	cout << "}" << endl;

	index_opt idx;
	search_interface *si = new mbes_opt(p, q, theta, sf, pf, sp, &ir, &idx);

	// reduction
	if (rd)
	{
		ir.clear_t();
		si->reduction_1hop(bg);
		si->reduction_2hop(bg);
		si->reduction_1hop(bg);
		cout << "after reduction:=>{";
		bg.print_info();
		cout << "}" << endl;
		reduction_time = ir.get_time();
	}

	// sort
	ir.clear_t();
	bg.reorder(sm);
	reorder_time = ir.get_time();

	// only execute reduction
	if (rd == 2)
	{
		// must provide 10 parameters
		cout << "time info:=>{" << reduction_time << "," << reorder_time << "}" << endl;
		bg.write(_argv[10]);
		return 0;
	}

	// running
	ir.clear_t();
	idx.init(bg.size());
	si->solve(bg);
	run_time = ir.get_time();

#ifdef index_on
	// simulate pre-order of prefix tree
	ir.clear_t();
	idx.remove_duplicate();
	idx_time = ir.get_time();
#endif

	cout << "result size:=>{" << idx.size() << "," << idx.attempts() << "}" << endl;
	cout << "time info:=>{" << reduction_time << "," << reorder_time << "," << run_time << "," << idx_time << "}" << endl;

	return 0;
}
