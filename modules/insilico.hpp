#ifndef INSILICO_HPP
#define INSILICO_HPP

#include <types/cellmodels/cellmodel.hpp>
#include <types/cipa_features.hpp>
#include <types/drug_block_input.hpp>
#include <types/cvar_input.hpp>
#include <types/parameter.hpp>
#include <vector>

// insilico.hpp
// Simulate in-silico drug simulation based on Dutta et. al. 2017.
// The end results are steady-state values in one of the paces
// among 1000 paces in the last 250 paces.

// Return true if EAD happened in the selected pace.

#include <vector>
using std::vector;

int insilico(double conc, const Drug_Row &hill, const Drug_Row &herg, const Parameter *p_param, Cipa_Features &p_features, short sample_id,
             const Cvar_Row *cvar = nullptr);

// load last states from the 1000 paces steady-state control simulation.
// If error occured, keep using the initial condition from the model.
void set_initial_condition(Cellmodel *p_cell, char *ic_file_name);

// functions that will be run at the end of cycle.
void end_of_cycle_funct(short *pace_count, double *inet_auc, Cellmodel *p_cell, const Parameter *p_param, Cipa_Features &p_features, Cipa_Features &temp_features,
                        FILE *fp_vmdebug);

bool get_dvmdt_repol_max(Cellmodel *p_cell, Cipa_Features &p_features, const Parameter *p_param, const double tcurr, const short pace_count);

#endif
