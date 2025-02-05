#ifndef INSILICO_CVAR_HPP
#define INSILICO_CVAR_HPP

#include <types/cellmodels/cellmodel.hpp>
#include <types/cipa_features.hpp>
#include <types/drug_block_input.hpp>
#include <types/cvar_input.hpp>
#include <types/parameter.hpp>
#include <vector>

// insilico_cvar.hpp
// Similar with insilico, but for conductance variability

#include <vector>
using std::vector;

void insilico_cvar(double conc, 
const Drug_Row& hill, const Drug_Row& herg, const Cvar_Row& cvar,
const Parameter *p_param, Cipa_Features &p_features, short sample_id);

// load last states from the 1000 paces steady-state control simulation.
// If error occured, keep using the initial condition from the model.
void cvar_set_initial_condition(Cellmodel *p_cell, char *ic_file_name);

// functions that will be run at the end of cycle.
void cvar_end_of_cycle_funct(short *pace_count, Cellmodel *p_cell, const Parameter *p_param, Cipa_Features &p_features, Cipa_Features &temp_features, FILE *fp_vmdebug);

bool cvar_get_dvmdt_repol_max(Cellmodel *p_cell, Cipa_Features &p_features, const Parameter *p_param, const double tcurr, const short pace_count) ;



#endif
