#include "drug_simulation_bench.hpp"

#include <mpi.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functions/inputoutput.hpp>
#include <iterator>
#include <map>
#include <string>
#include <types/cipa_features.hpp>
#include <types/cml_consts.hpp>
#include <types/cvar_input.hpp>
#include <types/drug_block_input.hpp>
#include <types/mpi_profile.hpp>
#include <vector>

#include "insilico.hpp"
#include "postprocessing.hpp"

using std::multimap;
using std::string;
using std::vector;

int drug_simulation_bench(const Parameter *p_param, multimap<double, string> &map_time_series, Cipa_Features &p_features, Drug_Block_Input &hill) {
  char buffer[900];
  int error_code;

  // parsing ic50 and hill constants file into variables
  get_data_from_file<Drug_Row, Drug_Block_Input>(p_param->hill_file, hill);
  error_code = check_drug_data_content(hill, p_param);
  if (error_code > 0) return 1;

  // parsing herg constants file into variables
  Drug_Block_Input herg;
  get_data_from_file<Drug_Row, Drug_Block_Input>(p_param->herg_file, herg);

  // parsing cvar constants file into variables
  Cvar_Input cvar;
  if (p_param->is_cvar == 1) {
    get_data_from_file<Cvar_Row, Cvar_Input>(p_param->cvar_file, cvar);
  }

#ifdef ORD_DYN_2017
  error_code = check_drug_data_content(herg, p_param);
  if (error_code > 0) return 1;
  if (herg.size() != hill.size()) {
    mpi_printf(0, "Data size not equal. Hill: %d hERG %d\n", hill.size(), herg.size());
    return 1;
  }
#endif

  // get concentration values from input
  // and make the directories of them.
  std::vector<double> concs;
  concs.push_back(0.);
  strncpy(buffer, p_param->concs, sizeof(buffer));
  char *token = strtok(buffer, ",");
  while (token != NULL) {  // begin data tokenizing
    concs.push_back(strtod(token, NULL));
    token = strtok(NULL, ",");
  }  // end data tokenizing

  if (MPI_Profile::rank == 0) create_concs_directories(concs, p_param->drug_name);
  MPI_Barrier(MPI_COMM_WORLD);

  // Progress tracking variables
  int total_samples = hill.size();
  int local_samples = 0;
  for (short i = MPI_Profile::rank; i < total_samples; i += MPI_Profile::size) {
    local_samples++;
  }
  int local_tasks = local_samples * (concs.size() - 1);
  if (MPI_Profile::rank == 0) local_tasks += 1; // Add control task for rank 0
  int tasks_completed = 0;
  double last_reported_progress = -1.0; // To track last reported progress percentage

  short sample_id = 0;
  int group_id = 0;
  bool is_ead;
  double inal_auc_control = 0.0;
  double ical_auc_control = 0.0;

#if defined POSTPROCESSING
  // Will be used in the postprocessing only.
  vector<double> vec_repol_states;
  char repol_states_file[255];
#endif

  if (MPI_Profile::rank == 0) {
#if defined POSTPROCESSING
      vec_repol_states.clear();
      mpi_printf(0, "Running control Postprocessing simulation started. Skipped the insilico process and read the steady state values....\n");
      snprintf(repol_states_file, sizeof(repol_states_file), "%s/%.2lf/%s_%.2lf_repol_states_smp%d_%s.csv", 
      p_param->repol_states_folder, 0.00, p_param->drug_name, 0.00, sample_id, p_param->user_name);
      error_code = set_initial_condition_postprocessing(repol_states_file, vec_repol_states);
      if(error_code != 0) return 1;
      p_features.repol_states.clear();
      p_features.repol_states.insert(p_features.repol_states.begin(), vec_repol_states.begin(), vec_repol_states.end());
#else
    mpi_printf(0, "Running control in-silico simulation...\n");
    insilico(0., hill[sample_id], herg[sample_id], p_param, p_features, sample_id);
#endif
    mpi_printf(0, "Running control postprocessing simulation...\n");
    postprocessing(0., inal_auc_control, ical_auc_control, hill[sample_id], herg[sample_id], p_param, p_features, sample_id, group_id);
    inal_auc_control = p_features.inal_auc;
    ical_auc_control = p_features.ical_auc;

    // Update and report progress for control simulation
    tasks_completed++;
    double local_progress = (double)tasks_completed / local_tasks * 100.0;
    mpi_printf(0, "Rank 0: %.0f%% complete (%d/%d tasks, control simulation done)\n", 
               local_progress, tasks_completed, local_tasks);

  }
  MPI_Barrier(MPI_COMM_WORLD);

  printf("Before INaL at rank %d: %lf\n", MPI_Profile::rank, inal_auc_control);
  printf("Before ICaL at rank %d: %lf\n", MPI_Profile::rank, ical_auc_control);
  MPI_Bcast(&inal_auc_control, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ical_auc_control, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  printf("After INaL at rank %d: %lf\n", MPI_Profile::rank, inal_auc_control);
  printf("After ICaL at rank %d: %lf\n", MPI_Profile::rank, ical_auc_control);


  for (sample_id = MPI_Profile::rank; sample_id < hill.size(); sample_id += MPI_Profile::size) {
    if (sample_id >= hill.size()) break;  // Ensure no out-of-bounds access

    for (short idx = 1; idx < concs.size(); idx++) {
      // postprocessing only read initial states at repolarization.
      // if the certain file cannot be read, skip to the next iteration.
#if defined POSTPROCESSING
        vec_repol_states.clear();
        mpi_printf(0, "Postprocessing simulation started. Skipped the insilico process and read the steady state values....\n");
        snprintf(repol_states_file, sizeof(repol_states_file), "%s/%.2lf/%s_%.2lf_repol_states_smp%d_%s.csv", 
          p_param->repol_states_folder, concs[idx], p_param->drug_name, concs[idx], sample_id, p_param->user_name);
        error_code = set_initial_condition_postprocessing(repol_states_file, vec_repol_states);
        if(error_code != 0) continue;
        p_features.repol_states.clear();
        p_features.repol_states.insert(p_features.repol_states.begin(), vec_repol_states.begin(), vec_repol_states.end());
#else
        mpi_printf(0, "insilico simulation started. After that, followed by postprocessing...\n");
        insilico(concs[idx], hill[sample_id], herg[sample_id], p_param, p_features, sample_id, &cvar[sample_id]);
#endif
        postprocessing(concs[idx], inal_auc_control, ical_auc_control, hill[sample_id], herg[sample_id], p_param, p_features, sample_id,
                            group_id, &cvar[sample_id]);

        // Update and report local progress
        tasks_completed++;
        double local_progress = (double)tasks_completed / local_tasks * 100.0;

        // Report progress at 10% increments to reduce output overhead
        if (local_progress >= last_reported_progress + 10.0 || tasks_completed == local_tasks) {
          last_reported_progress = local_progress - fmod(local_progress, 10.0); // Round down to nearest 10%
          printf("Rank %d: %.0f%% complete (%d/%d tasks)\n", 
                     MPI_Profile::rank, local_progress, tasks_completed, local_tasks);
        }

        // Aggregate progress across all ranks
        double global_progress = 0.0;
        double local_progress_fraction = (double)tasks_completed / local_tasks;
        MPI_Reduce(&local_progress_fraction, &global_progress, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        // Rank 0 computes and displays total progress
        if (MPI_Profile::rank == 0) {
          global_progress = (global_progress / MPI_Profile::size) * 100.0;
          mpi_printf(0, "DrugSimulation Global progress: %.0f%%\n", global_progress);
        }
    }
    
  }
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}


#if defined POSTPROCESSING
int set_initial_condition_postprocessing(const char *repol_states_file, vector<double> &vec_repol_states) {
  char buffer[50];
  FILE *fp_repol_states;

  fp_repol_states = fopen(repol_states_file, "r");
  if (fp_repol_states == NULL) {
    mpi_printf(cml::commons::MASTER_NODE, "File %s not found! Make sure the name is correct!\n", repol_states_file);
    return 1;
  }
  mpi_printf(cml::commons::MASTER_NODE, "Using states from repolarization phase from previous simulations\n");
  while (fgets(buffer, sizeof(buffer), fp_repol_states) != NULL) {
    vec_repol_states.push_back(strtod(buffer, NULL));
  }
  fclose(fp_repol_states); 

  return 0;
}
#endif
