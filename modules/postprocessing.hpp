#ifndef POSTPROCESSING_HPP
#define POSTPROCESSING_HPP

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

void postprocessing(double conc, double inal_auc_control, double ical_auc_control, const Drug_Row &hill, const Drug_Row &herg,
                    const Parameter *p_param, Cipa_Features &p_features, short sample_id, short group_id, const Cvar_Row *cvar = nullptr);

// load last states from the 1000 paces steady-state control simulation.
// If error occured, keep using the initial condition from the model.
short set_initial_condition_postprocessing(Cellmodel *p_cell, const char *ic_file_name);

// functions related to the features retrieval
void get_vm_features_postprocessing(Cellmodel *p_cell, Cipa_Features &p_features, const double tcurr);
void get_ca_features_postprocessing(Cellmodel *p_cell, Cipa_Features &p_features, const double tcurr);
void get_duration_postprocessing(Cipa_Features &p_features);
void collect_features(Cipa_Features &p_features, const Parameter *p_param, Cellmodel *p_cell, double conc, double inet_auc, double inet_apd_auc,
                      double inal_auc, double ical_auc, double inal_auc_control, double ical_auc_control, char *features_file_name, short sample_id,
                      short group_id);
#endif
