#include "insilico.hpp"

#ifdef TOMEK_2019
#include <types/cellmodels/Tomek_model.hpp>
#elif defined(TOMEK_DYNCL_2020)
#include <types/cellmodels/Tomek_dynCl.hpp>
#elif defined(ORD_DYN_2017)
#include <types/cellmodels/ohara_rudy_cipa_v1_2017.hpp>
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

void insilico(double conc, const Drug_Row &hill, const Drug_Row &herg, const Parameter *p_param, Cipa_Features &p_features, short sample_id,
              const Cvar_Row *cvar) {
  // simulation parameters.
  // Assigned to the constant variables for
  // the sake of simplicity.
  const double bcl = p_param->bcl;
  const short pace_max = p_param->pace_max;
  const short celltype = p_param->celltype;
  const char *user_name = p_param->user_name;
  const double dt_min = p_param->dt_min;
  const double dt_max = p_param->dt_max;
  const double dtw = p_param->dtw;
  const double stim_dur = p_param->stim_dur;
  const double stim_amp_scale = p_param->stim_amp_scale;
  const short solver_type = p_param->solver_type;
  const short is_cvar = p_param->is_cvar;
  const char *drug_name = p_param->drug_name;

  // this is the cellmodel initialization part
  Cellmodel *p_cell;
#if defined ORD_DYN_2017
  mpi_printf(cml::commons::MASTER_NODE, "Using ORd-dyn 2017 cell model\n");
  p_cell = new ohara_rudy_cipa_v1_2017();
  if (is_cvar && cvar) {
    p_cell->initConsts(static_cast<double>(celltype), conc, hill.data, herg.data, cvar->data);
  } else {
    p_cell->initConsts(static_cast<double>(celltype), conc, hill.data, herg.data);
  }
#elif defined TOMEK_2019
  mpi_printf(cml::commons::MASTER_NODE, "Using Tomek 2019 cell model\n");
  p_cell = new Tomek_model();
  if (is_cvar && cvar) {
    p_cell->initConsts(static_cast<double>(celltype), conc, hill.data, cvar->data);
  } else {
    p_cell->initConsts(static_cast<double>(celltype), conc, hill.data);
  }
#elif defined TOMEK_DYNCL_2020
  mpi_printf(cml::commons::MASTER_NODE, "Using Tomek_dynCl 2020 cell model\n");
  p_cell = new Tomek_dynCl();
  if (is_cvar && cvar) {
    p_cell->initConsts(static_cast<double>(celltype), conc, hill.data, cvar->data);
  } else {
    p_cell->initConsts(static_cast<double>(celltype), conc, hill.data);
  }
#else
  mpi_printf(cml::commons::MASTER_NODE, "Using ORd2011 cell model\n");
  p_cell = new Ohara_Rudy_2011();
  if (is_cvar && cvar) {
    p_cell->initConsts(static_cast<double>(celltype), conc, hill.data, true, cvar->data);
  } else {
    p_cell->initConsts(static_cast<double>(celltype), conc, hill.data, true);
  }
#endif

  p_cell->CONSTANTS[BCL] = bcl;
  p_cell->CONSTANTS[duration] = stim_dur;
  p_cell->CONSTANTS[amp] *= stim_amp_scale;

  // variables for I/O
  char buffer[900];
  FILE *fp_vmdebug, *fp_last_states;

  snprintf(buffer, sizeof(buffer), "%s/%s/%.2lf/%s_%.2lf_vmdebug_smp%d_%s.csv", cml::commons::RESULT_FOLDER, drug_name, conc, drug_name, conc, sample_id, user_name);
  fp_vmdebug = fopen(buffer, "w");
  fprintf(fp_vmdebug, "%s,%s,%s,%s,%s,%s,%s,%s\n", "Pace", "T_Peak", "Vmpeak", "Vmvalley", "Vm_repol30", "Vm_repol50", "Vm_repol90",
          "dVmdt_repol_max");
  snprintf(buffer, sizeof(buffer), "%s/%s/%.2lf/%s_%.2lf_last_states_smp%d_%s.csv", cml::commons::RESULT_FOLDER, drug_name, conc, drug_name, conc, sample_id, user_name);
  fp_last_states = fopen(buffer, "w");

  FILE *fp_last_ten_paces;
  snprintf(buffer, sizeof(buffer), "%s/%s/%.2lf/%s_%.2lf_last_10_paces_smp%d_%s.csv", cml::commons::RESULT_FOLDER, drug_name, conc, drug_name, conc, sample_id, user_name);
  fp_last_ten_paces = fopen(buffer, "w");
  fprintf(fp_last_ten_paces, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", "Time", "Vm", "dVm/dt", "timestep(dt)", "Cai(x1.000.000)(milliM->nanoM)",
          "INa(x1.000)(microA->nanoA)", "INaL(x1.000)(microA->nanoA)", "ICaL(x1.000)(microA->nanoA)", "Ito(x1.000)(microA->nanoA)",
          "IKr(x1.000)(microA->nanoA)", "IKs(x1.000)(microA->nanoA)", "IK1(x1.000)(microA->nanoA)");

  snprintf(buffer, sizeof(buffer), "last_states_1000paces_%s.dat", p_param->cell_name);
  mpi_printf(cml::commons::MASTER_NODE, "Last steady-state file: %s\n", buffer);
  // replace the initial condition
  // with the last state value from
  // steady-state 1000 paces control simulation.
  // set_initial_condition(p_cell, buffer);

  // time variables.
  // some of these values are taken from the supplementary materials of ORd2011
  // page 19.
  double dt = p_param->dt;
  double dt_set;
  double tcurr = 0.;
  double tmax = pace_max * bcl;
  double time_point = 25.0;
  double tprint = 0.;
  short pace_count = 0;
  long icount = 0;
  long imax = (long)((pace_max * bcl)/dt);
  const long print_freq = (int)(1./dt) * dtw;
  double next_print_time = 0.0;

  // CVode solver.
  CVodeSolverData *p_cvode;
  int cvode_retval;
  if (solver_type == 0) {
    mpi_printf(cml::commons::MASTER_NODE, "Using Euler Solver.\n");
  } else if (solver_type == 1) {
    mpi_printf(cml::commons::MASTER_NODE, "Using CVode Solver.\n");
    p_cvode = new CVodeSolverData();
    init_cvode(p_cvode, p_cell, tcurr);
    set_dt_cvode(p_cvode, tcurr, time_point, bcl, dt_min, dt_max, &dt);
  }

  Cipa_Features temp_features;
  temp_features.init(p_cell->STATES[V], p_cell->STATES[cai] * cml::math::MILLI_TO_NANO);
  p_features.init(p_cell->STATES[V], p_cell->STATES[cai] * cml::math::MILLI_TO_NANO);

  // begin computation loop
  while (icount < imax) {
    // compute and solving part
    //
    // Different solver_type has different
    // method of calling computeRates().
    if (solver_type == 0) {
      p_cell->computeRates(tcurr, p_cell->CONSTANTS, p_cell->RATES, p_cell->STATES, p_cell->ALGEBRAIC);
      solveEuler(dt, p_cell);
      // increase the time based on the dt.
      tcurr += dt;
    } else if (solver_type == 1) {
      cvode_retval = solve_cvode(p_cvode, p_cell, tcurr + dt, &tcurr);
      if (cvode_retval != 0) {
        printf("CVode calculation error happened at sample %hd concentration %.0lf!!\n", sample_id, conc);
        return;
      }
      set_dt_cvode(p_cvode, tcurr, time_point, bcl, dt_min, dt_max, &dt);
    }

    // execute this code when entering new cycle.
    if (floor(tcurr / bcl) != pace_count) {
      end_of_cycle_funct(&pace_count, p_cell, p_param, p_features, temp_features, fp_vmdebug);
      // mpi_printf(0,"Entering pace %hd at %lf msec.\n", pace_count, tcurr);
    }

    // begin the last 250 paces operations.
    if (pace_count >= pace_max - 250) {
      get_dvmdt_repol_max(p_cell, temp_features, p_param, tcurr, pace_count);
    }  // end of the last 250 paces operations.

    // We use 1000. as the multiplier
    // assuming that the minimum dt is 0.005.
    // Because of this,
    // we will have dtw universal for any solver.
    if (tcurr >= next_print_time) {
    //if (icount% print_freq == 0) {
      if (tcurr > bcl * (pace_max - 10)) {
        tprint = tcurr - (bcl * (pace_max - 10));
        snprintf(buffer, sizeof(buffer), "%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", p_cell->STATES[V], p_cell->RATES[V],
                 dt, p_cell->STATES[cai] * cml::math::MILLI_TO_NANO, p_cell->ALGEBRAIC[INa] * cml::math::MICRO_TO_NANO,
                 p_cell->ALGEBRAIC[INaL] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[ICaL] * cml::math::MICRO_TO_NANO,
                 p_cell->ALGEBRAIC[Ito] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[IKr] * cml::math::MICRO_TO_NANO,
                 p_cell->ALGEBRAIC[IKs] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[IK1] * cml::math::MICRO_TO_NANO);
        fprintf(fp_last_ten_paces, "%.0lf,%s", round(tprint), buffer);
      }
      next_print_time += dtw;
    }
    icount++;

  }  // end computation loop

  // write the last states into a file for further usage.
  for (short idx = 0; idx < p_cell->states_size; idx++) {
    fprintf(fp_last_states, "%.16lf\n", p_features.last_states[idx]);
  }

  fprintf(fp_vmdebug, "Selected final pace: %hd\n", p_features.pace_target);
  fprintf(fp_vmdebug, "Features saved: \n%s,%s,%s,%s,%s,%s,%s,%s\n", "Pace", "T_Peak", "Vmpeak", "Vmvalley", "Vm_repol30", "Vm_repol50", "Vm_repol90",
          "dVmdt_repol_max");
  fprintf(fp_vmdebug, "%hd,%.4lf,%.4lf,%.4lf,%.4lf,%.8lf,%.8lf,%.8lf\n", p_features.pace_target, p_features.time_vm_peak, p_features.vm_peak,
          p_features.vm_valley, p_features.vm_amp30, p_features.vm_amp50, p_features.vm_amp90, p_features.dvmdt_repol_max);
  fprintf(fp_vmdebug, "Last states:\n");
  for (short idx = 0; idx < p_cell->states_size; idx++) {
    fprintf(fp_vmdebug, "%.16lf\n", p_features.last_states[idx]);
  }

  delete p_cvode;
  fclose(fp_last_ten_paces);
  fclose(fp_last_states);
  fclose(fp_vmdebug);
  
}

void end_of_cycle_funct(short *pace_count, Cellmodel *p_cell, const Parameter *p_param, Cipa_Features &p_features, Cipa_Features &temp_features,
                        FILE *fp_vmdebug) {
  if (*pace_count >= p_param->pace_max - 250) {
    fprintf(fp_vmdebug, "%hd,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", *pace_count, temp_features.time_vm_peak, temp_features.vm_peak,
            temp_features.vm_valley, temp_features.vm_amp30, temp_features.vm_amp50, temp_features.vm_amp90, temp_features.dvmdt_repol_max);
    fflush(fp_vmdebug);
    temp_features.last_states.insert(temp_features.last_states.begin(), p_cell->STATES, p_cell->STATES + p_cell->states_size);
    /*
        mpi_printf(cml::commons::MASTER_NODE, "Last states of pace %d:\n", *pace_count);
        for(short idx = 0; idx < p_cell->states_size; idx++){
          mpi_printf(cml::commons::MASTER_NODE, "%lf,", temp_features.last_states[idx]);
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
    mpi_printf(cml::commons::MASTER_NODE, "No initial file found! Keep using initial condition from the cell model!\n");
  }
}

bool get_dvmdt_repol_max(Cellmodel *p_cell, Cipa_Features &p_features, const Parameter *p_param, const double tcurr, const short pace_count) {
  bool is_eligible_AP = false;
  static const double TIME_NORMALIZER = (p_param->bcl * (p_param->pace_max - 250));

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
        p_features.time_vm_peak = fmod((tcurr - TIME_NORMALIZER), p_param->bcl);
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

