#include "insilico.hpp"

#if defined CIPAORDV1_0
#include <types/cellmodels/ohara_rudy_cipa_v1_2017.hpp>
#elif defined TOR_ORD
#include <types/cellmodels/Tomek_model.hpp>
#elif defined TOR_ORD_DYNCL
#include <types/cellmodels/Tomek_dynCl.hpp>
#else
#include <types/cellmodels/Ohara_Rudy_2011.hpp>
#endif

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functions/helper_cvode.hpp>
#include <functions/helper_math.hpp>
#include <functions/inputoutput.hpp>
#include <map>
#include <types/cipa_features.hpp>
#include <types/cml_consts.hpp>
#include <vector>

using std::vector;

int insilico(double conc, const Drug_Row &hill, const Drug_Row &herg, const Parameter *p_param, Cipa_Features &p_features, short sample_id,
              const Cvar_Row *cvar) {
  // simulation parameters.
  // Assigned to the constant variables for
  // the sake of simplicity.
  const double cycle_length = p_param->cycle_length;
  const short number_pacing = p_param->number_pacing;
  const short number_pacing_write = p_param->number_pacing_write;
  const char *user_name = p_param->user_name;
  const double time_step_min = p_param->time_step_min;
  const double time_step_max = p_param->time_step_max;
  const double writing_step = p_param->writing_step;
  const double stimulus_duration = p_param->stimulus_duration;
  const double stimulus_amplitude_scale = p_param->stimulus_amplitude_scale;
  const char *solver_type = p_param->solver_type;
  const short is_cvar = p_param->is_cvar;
  const char *drug_name = p_param->drug_name;

  if(time_step_min > writing_step){
    mpi_printf(cml::commons::MASTER_NODE,"%s\n%s\n",
    "WARNING!!! The writing_step values is smaller than the timestep!",
    "Simulation still run, but the time series will use time_step_min as writing step.");
  }

  // this is the cellmodel initialization part
  Cellmodel *p_cell;
  short cell_type;
  const char *cell_model = p_param->cell_model;
  if( strstr(cell_model,"endo") != NULL ) cell_type = 0;
  else if( strstr(cell_model,"epi") != NULL ) cell_type = 1;
  else if( strstr(cell_model,"myo") != NULL ) cell_type = 2;
  mpi_printf(cml::commons::MASTER_NODE,"Using %s cell model with cell_type=%hd\n", cell_model, cell_type);
#if defined CIPAORDV1_0
  p_cell = new ohara_rudy_cipa_v1_2017();
  if (is_cvar && cvar) {
    p_cell->initConsts(static_cast<double>(cell_type), conc, hill.data, herg.data, cvar->data);
  } else {
    p_cell->initConsts(static_cast<double>(cell_type), conc, hill.data, herg.data);
  }
#elif defined TOR_ORD
  p_cell = new Tomek_model();
  if (is_cvar && cvar) {
    p_cell->initConsts(static_cast<double>(cell_type), conc, hill.data, cvar->data);
  } else {
    p_cell->initConsts(static_cast<double>(cell_type), conc, hill.data);
  }
#elif defined TOR_ORD_DYNCL
  p_cell = new Tomek_dynCl();
  if (is_cvar && cvar) {
    p_cell->initConsts(static_cast<double>(cell_type), conc, hill.data, cvar->data);
  } else {
    p_cell->initConsts(static_cast<double>(cell_type), conc, hill.data);
  }
#else
  p_cell = new Ohara_Rudy_2011();
  if (is_cvar && cvar) {
    p_cell->initConsts(static_cast<double>(cell_type), conc, hill.data, true, cvar->data);
  } else {
    p_cell->initConsts(static_cast<double>(cell_type), conc, hill.data, true);
  }
#endif

  p_cell->CONSTANTS[BCL] = cycle_length;
  p_cell->CONSTANTS[duration] = stimulus_duration;
  p_cell->CONSTANTS[amp] *= stimulus_amplitude_scale;

  // variables for I/O
  char buffer[900];
  FILE *fp_vmdebug, *fp_initial_values, *fp_last_paces;

  snprintf(buffer, sizeof(buffer), "%s/%.2lf/%s_%.2lf_vmdebug_smp%d.csv", 
          cml::commons::RESULT_FOLDER, conc, drug_name, conc, sample_id);
  fp_vmdebug = fopen(buffer, "w");
  if(fp_vmdebug == NULL){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr, "Cannot create file %s. Make sure the directory is existed!!!\n",buffer);
    return 1;
  }
  fprintf(fp_vmdebug, "%s,%s,%s,%s,%s,%s,%s,%s\n", "Pace", "T_Peak", "Vmpeak", "Vmvalley", "Vm_repol30", "Vm_repol50", "Vm_repol90",
          "dVmdt_repol_max");

  snprintf(buffer, sizeof(buffer), "%s/%.2lf/%s_%.2lf_initial_values_smp%d.csv", 
          cml::commons::RESULT_FOLDER, conc, drug_name, conc, sample_id);
  fp_initial_values = fopen(buffer, "w");
  if(fp_initial_values == NULL){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr, "Cannot create file %s. Make sure the directory is existed!!!\n",buffer);
    return 1;
  }

  snprintf(buffer, sizeof(buffer), "%s/%.2lf/%s_%.2lf_last_paces_smp%hd.csv", 
          cml::commons::RESULT_FOLDER, conc, drug_name, conc, sample_id);
  fp_last_paces = fopen(buffer, "w");
  if(fp_last_paces == NULL){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr, "Cannot create file %s. Make sure the directory is existed!!!\n",buffer);
    return 1;
  }
  fprintf(fp_last_paces, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", 
          "Time(ms)", "Last_Vm(mV)", "Last_dVm/dt(mV/ms)", "Last_Cai(nM)",
          "Last_INa(nA)", "Last_INaL(nA)", "Last_ICaL(nA)", "Last_Ito(nA)",
          "Last_IKr(nA)", "Last_IKs(nA)", "Last_IK1(nA)");

  //snprintf(buffer, sizeof(buffer), "./%s/last_states_%hdpaces_%s.dat", "initial_values_files", number_pacing, p_param->cell_model);
  //mpi_printf(cml::commons::MASTER_NODE, "Last steady-state file: %s\n", buffer);
  // replace the initial condition
  // with the last state value from
  // steady-state 1000 paces control simulation.
  // TODO: implemented later on.
  //set_initial_condition(p_cell, buffer);

  // time variables.
  // some of these values are taken from the supplementary materials of ORd2011
  // page 19.
  double time_step = time_step_min;
  double tcurr = 0.;
  double tmax = number_pacing * cycle_length;
  double time_point = 25.0;
  double tprint = 0.;
  short pace_count = 0;

  // CVode solver.
  CVodeSolverData *p_cvode;
  int cvode_retval;
  if (strncasecmp(solver_type, "Euler", sizeof(solver_type)) == 0) {
    mpi_printf(cml::commons::MASTER_NODE, "Using Euler Solver.\n");
  } else if (strncasecmp(solver_type, "CVode", sizeof(solver_type)) == 0) {
    mpi_printf(cml::commons::MASTER_NODE, "Using CVode Solver.\n");
    p_cvode = new CVodeSolverData();
    init_cvode(p_cvode, p_cell, tcurr);
    set_dt_cvode(p_cvode, tcurr, time_point, cycle_length, time_step_min, time_step_max, &time_step);
  }
  else{
    mpi_fprintf(cml::commons::MASTER_NODE, stderr, "Solver type %s is undefined! Please choose the available solver types from the manual!\n", solver_type);
    return 1;
  }


  Cipa_Features temp_features;
  temp_features.init(p_cell->STATES[V], p_cell->STATES[cai] * cml::math::MILLI_TO_NANO);
  p_features.init(p_cell->STATES[V], p_cell->STATES[cai] * cml::math::MILLI_TO_NANO);

  // begin computation loop
  while (tcurr < tmax) {
    // compute and solving part
    //
    // Different solver_type has different
    // method of calling computeRates().
    if (strncasecmp(solver_type, "Euler", sizeof(solver_type)) == 0) {
      p_cell->computeRates(tcurr, p_cell->CONSTANTS, p_cell->RATES, p_cell->STATES, p_cell->ALGEBRAIC);
      solveEuler(time_step, p_cell);
      // increase the time based on the dt.
      tcurr += time_step;
    } else if (strncasecmp(solver_type, "CVode", sizeof(solver_type)) == 0) {
      cvode_retval = solve_cvode(p_cvode, p_cell, tcurr + time_step, &tcurr);
      if (cvode_retval != 0) {
        printf("CVode calculation error happened at sample %hd concentration %.0lf!!\n", sample_id, conc);
        return cvode_retval;
      }
      set_dt_cvode(p_cvode, tcurr, time_point, cycle_length, time_step_min, time_step_max, &time_step);
    }

    // execute this code when entering new cycle.
    if (floor(tcurr / cycle_length) != pace_count) {
      end_of_cycle_funct(&pace_count, p_cell, p_param, p_features, temp_features, fp_vmdebug);
      // mpi_printf(0,"Entering pace %hd at %lf msec.\n", pace_count, tcurr);
    }

    // begin the last 250 paces operations.
    if (pace_count >= number_pacing - 250) {
      get_dvmdt_repol_max(p_cell, temp_features, p_param, tcurr, pace_count);
    }  // end of the last 250 paces operations.

    // Output part.
    if (tcurr > cycle_length * (number_pacing - number_pacing_write)) {
      tprint = tcurr - (cycle_length * (number_pacing - number_pacing_write));
      snprintf(buffer, sizeof(buffer), "%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", p_cell->STATES[V], p_cell->RATES[V],
               p_cell->STATES[cai] * cml::math::MILLI_TO_NANO, p_cell->ALGEBRAIC[INa] * cml::math::MICRO_TO_NANO,
               p_cell->ALGEBRAIC[INaL] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[ICaL] * cml::math::MICRO_TO_NANO,
               p_cell->ALGEBRAIC[Ito] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[IKr] * cml::math::MICRO_TO_NANO,
               p_cell->ALGEBRAIC[IKs] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[IK1] * cml::math::MICRO_TO_NANO);
      fprintf(fp_last_paces, "%.4lf,%s", tprint, buffer);
    }
    tprint += writing_step;

  }  // end computation loop

  // write the last states into a file for further usage.
  for (short idx = 0; idx < p_cell->states_size; idx++) {
    fprintf(fp_initial_values, "%.16lf\n", p_features.initial_values[idx]);
  }

  fprintf(fp_vmdebug, "Selected final pace: %hd\n", p_features.pace_target);
  fprintf(fp_vmdebug, "Features saved: \n%s,%s,%s,%s,%s,%s,%s,%s\n", "Pace", "T_Peak", "Vmpeak", "Vmvalley", "Vm_repol30", "Vm_repol50", "Vm_repol90",
          "dVmdt_repol_max");
  fprintf(fp_vmdebug, "%hd,%.4lf,%.4lf,%.4lf,%.4lf,%.8lf,%.8lf,%.8lf\n", p_features.pace_target, p_features.time_vm_peak, p_features.vm_peak,
          p_features.vm_valley, p_features.vm_amp30, p_features.vm_amp50, p_features.vm_amp90, p_features.dvmdt_repol_max);
  fprintf(fp_vmdebug, "States during the highest dVM/dt between 30%% and 90%% repolarization phase:\n");
  for (short idx = 0; idx < p_cell->states_size; idx++) {
    fprintf(fp_vmdebug, "%.8lf\n", p_features.initial_values[idx]);
  }

  if( strncasecmp(solver_type, "CVode", sizeof(solver_type)) == 0 ){
    clean_cvode(p_cvode);
    delete p_cvode;
  }
  fclose(fp_last_paces);
  fclose(fp_initial_values);
  fclose(fp_vmdebug);
  
  return 0;
}

void end_of_cycle_funct(short *pace_count, Cellmodel *p_cell, const Parameter *p_param, Cipa_Features &p_features, Cipa_Features &temp_features,
                        FILE *fp_vmdebug) {
  if (*pace_count >= p_param->number_pacing - 250) {
    fprintf(fp_vmdebug, "%hd,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", *pace_count, temp_features.time_vm_peak, temp_features.vm_peak,
            temp_features.vm_valley, temp_features.vm_amp30, temp_features.vm_amp50, temp_features.vm_amp90, temp_features.dvmdt_repol_max);
    fflush(fp_vmdebug);
    temp_features.initial_values.insert(temp_features.initial_values.begin(), p_cell->STATES, p_cell->STATES + p_cell->states_size);
    /*
        mpi_printf(cml::commons::MASTER_NODE, "Last states of pace %d:\n", *pace_count);
        for(short idx = 0; idx < p_cell->states_size; idx++){
          mpi_printf(cml::commons::MASTER_NODE, "%lf,", temp_features.initial_values[idx]);
        }
        mpi_printf(cml::commons::MASTER_NODE, "\n");
    */
    // mpi_printf(cml::commons::MASTER_NODE,"Temp dvmdt_repol_max: %lf ---- current dvmdt_repol_max: %lf\n", temp_features.dvmdt_repol_max,
    // p_features.dvmdt_repol_max);
    if (temp_features.dvmdt_repol_max > p_features.dvmdt_repol_max) {
      temp_features.pace_target = *pace_count;
      p_features = temp_features;
      // mpi_printf(cml::commons::MASTER_NODE,"Selected pace: %hd\n", p_features.pace_target);
    }
  }

  *pace_count += 1;
  temp_features.init(p_cell->STATES[V], p_cell->STATES[cai] * cml::math::MILLI_TO_NANO);
}

void set_initial_condition(Cellmodel *p_cell, char *ic_file_name) {
  char buffer[50];
  FILE *fp_states;
  short idx;
  mpi_printf(cml::commons::MASTER_NODE, "STATES before:\n");
  for (idx = 0; idx < 10; idx++) {
    mpi_printf(cml::commons::MASTER_NODE, "%lf ", p_cell->STATES[idx]);
  }
  mpi_printf(cml::commons::MASTER_NODE, "\n");
  fp_states = fopen(ic_file_name, "r");
  if (fp_states != NULL) {
    mpi_printf(cml::commons::MASTER_NODE, "Using initial condition from 1000 paces steady-state!\n");
    idx = 0;
    while (fgets(buffer, sizeof(buffer), fp_states) != NULL) {
      p_cell->STATES[idx++] = strtod(buffer, NULL);
    }
    mpi_printf(cml::commons::MASTER_NODE, "STATES after:\n");
    for (idx = 0; idx < 10; idx++) {
      mpi_printf(cml::commons::MASTER_NODE, "%lf ", p_cell->STATES[idx]);
    }
    mpi_printf(cml::commons::MASTER_NODE, "\n");
    fclose(fp_states);
  } else {
    mpi_printf(cml::commons::MASTER_NODE, "File %s\n cannot be found! Keep using initial condition from the cell model!\n", buffer);
  }
}

bool get_dvmdt_repol_max(Cellmodel *p_cell, Cipa_Features &p_features, const Parameter *p_param, const double tcurr, const short pace_count) {
  bool is_eligible_AP = false;
  static const double TIME_NORMALIZER = (p_param->cycle_length * (p_param->number_pacing - 250));

  // Find peak vm around 10 msecs and  30 msecs after stimulation
  if (tcurr > ((p_cell->CONSTANTS[BCL] * pace_count) + (p_cell->CONSTANTS[stim_start] + 10)) &&
      tcurr < ((p_cell->CONSTANTS[BCL] * pace_count) + (p_cell->CONSTANTS[stim_start] + 30))) {
    if (p_cell->STATES[V] > p_features.vm_peak) {
      p_features.vm_peak = p_cell->STATES[V];
      if (p_features.vm_peak > 0) {
        is_eligible_AP = true;
        p_features.vm_amp30 = p_features.vm_peak - (0.3 * (p_features.vm_peak - p_features.vm_valley));
        p_features.vm_amp50 = p_features.vm_peak - (0.5 * (p_features.vm_peak - p_features.vm_valley));
        p_features.vm_amp90 = p_features.vm_peak - (0.9 * (p_features.vm_peak - p_features.vm_valley));
        p_features.time_vm_peak = fmod((tcurr - TIME_NORMALIZER), p_param->cycle_length);
      }
    }
  }

  if (tcurr > ((p_cell->CONSTANTS[BCL] * pace_count) + (p_cell->CONSTANTS[stim_start] + 30))) {
    if (p_cell->RATES[V] > p_features.dvmdt_repol_max && p_features.vm_amp30 >= p_cell->STATES[V] && p_cell->STATES[V] >= p_features.vm_amp90) {
      p_features.dvmdt_repol_max = p_cell->RATES[V];
    }
  }

  return true;
}

