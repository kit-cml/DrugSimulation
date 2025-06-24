#ifndef DRUG_SIMULATION_BENCH_HPP
#define DRUG_SIMULATION_BENCH_HPP

#include <types/drug_block_input.hpp>
#include <types/parameter.hpp>
#include <types/cipa_features.hpp>

#include <map>
#include <string>

using std::multimap;
using std::string;

int drug_simulation_bench(const Parameter *p_param, multimap<double, string> &map_time_series, Cipa_Features &p_features, Drug_Block_Input &hill);
#if defined POSTPROCESSING
int set_initial_condition_postprocessing(const char *repol_states_file, vector<double> &vec_repol_states);
#endif

#endif
