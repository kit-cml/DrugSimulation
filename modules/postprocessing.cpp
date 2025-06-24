#include "postprocessing.hpp"

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
#include <fstream>
#include <functions/helper_cvode.hpp>
#include <functions/helper_math.hpp>
#include <functions/inputoutput.hpp>
#include <map>
#include <types/cipa_features.hpp>
#include <types/cml_consts.hpp>
#include <types/mpi_profile.hpp>
#include <vector>

using std::copy;
using std::vector;

void postprocessing(double conc, double inal_auc_control, double ical_auc_control, const Drug_Row &hill, const Drug_Row &herg,
                    const Parameter *p_param, Cipa_Features &p_features, short sample_id, short group_id, const Cvar_Row *cvar) {
  bool is_ead = false;
  // simulation parameters.
  // Assigned to the constant variables for
  // the sake of simplicity.
  const double bcl = p_param->bcl;
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
  FILE *fp_time_series;

  snprintf(buffer, sizeof(buffer), "%s/%s/%.2lf/%s_%.2lf_repol_states_smp%d_%s.csv", cml::commons::RESULT_FOLDER, drug_name, conc, drug_name, conc, sample_id, user_name);
  mpi_printf(cml::commons::MASTER_NODE, "Last steady-state file: %s\n", buffer);
  // replace the initial condition
  // with the last state value from
  // the in-silico simulation
  mpi_printf(cml::commons::MASTER_NODE, "STATES before:\n");
  for (short idx = 0; idx < 20; idx++) {
    mpi_printf(cml::commons::MASTER_NODE, "%lf ", p_cell->STATES[idx]);
  }
  mpi_printf(cml::commons::MASTER_NODE, "\nSTATES from repol_states vector:\n");
  for (short idx = 0; idx < 20; idx++) {
    mpi_printf(cml::commons::MASTER_NODE, "%lf ", p_features.repol_states[idx]);
  }
  mpi_printf(cml::commons::MASTER_NODE, "\nUsing repol state from the in-silico simulation.\n");

  copy(p_features.repol_states.begin(), p_features.repol_states.end(), p_cell->STATES);
  // set_initial_condition_postprocessing(p_cell, buffer);

  mpi_printf(cml::commons::MASTER_NODE, "STATES after:\n");
  for (short idx = 0; idx < 10; idx++) {
    mpi_printf(cml::commons::MASTER_NODE, "%lf ", p_cell->STATES[idx]);
  }
  mpi_printf(cml::commons::MASTER_NODE, "\n");

  snprintf(buffer, sizeof(buffer), "%s/%s/%.2lf/%s_%.2lf_time_series_smp%d_%s.csv", cml::commons::RESULT_FOLDER, drug_name, conc, drug_name, conc, sample_id, user_name);
  fp_time_series = fopen(buffer, "w");
  fprintf(fp_time_series, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", "Time(msec)", "Vm(mVolt)", "dVm/dt(mVolt/msec)", "timestep(dt)",
          "Cai(x1.000.000)(milliM->nanoM)", "INa(x1.000)(microA->nanoA)", "INaL(x1.000)(microA->nanoA)", "ICaL(x1.000)(microA->nanoA)",
          "Ito(x1.000)(microA->nanoA)", "IKr(x1.000)(microA->nanoA)", "IKs(x1.000)(microA->nanoA)", "IK1(x1.000)(microA->nanoA)", "Inet(microA)",
          "Inet_APD(microA)");

  // time variables.
  // some of these values are taken from the supplementary materials of ORd2011
  // page 19.
  double dt = p_param->dt;
  double dt_set;
  double tcurr = 0.;
  double time_point = 25.;
  double tmax = bcl;
  short pace_count = 0;
  long icount = 0;
  long imax = (long)((bcl)/dt);
  const long print_freq = (long)(1./dt) * dtw;
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

  double inal_auc = 0.;
  double ical_auc = 0.;
  double inet = 0.;
  double inet_apd = 0.;
  double inet_auc = 0.;
  double inet_apd_auc = 0.;
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

    // calculating the AUC
    inet = p_cell->ALGEBRAIC[INaL] + p_cell->ALGEBRAIC[ICaL] + p_cell->ALGEBRAIC[Ito] + p_cell->ALGEBRAIC[IKr] + p_cell->ALGEBRAIC[IKs] +
           p_cell->ALGEBRAIC[IK1];
    inet_auc += inet * dt;
    if (p_cell->STATES[V] > p_features.vm_amp90) {
      inet_apd = p_cell->ALGEBRAIC[INaL] + p_cell->ALGEBRAIC[ICaL] + p_cell->ALGEBRAIC[Ito] + p_cell->ALGEBRAIC[IKr] + p_cell->ALGEBRAIC[IKs] +
                 p_cell->ALGEBRAIC[IK1];
      inet_apd_auc += inet_apd * dt;
    }
    inal_auc += p_cell->ALGEBRAIC[INaL] * dt;
    ical_auc += p_cell->ALGEBRAIC[ICaL] * dt;

    // finding other vm and cai features
    get_vm_features_postprocessing(p_cell, p_features, tcurr);
    get_ca_features_postprocessing(p_cell, p_features, tcurr);

    // write the result to the graph
    if (tcurr >= next_print_time) {
    //if (icount% print_freq == 0) {
      // mpi_printf(0,"Writing at %lf msec.\n", tcurr);
      snprintf(buffer, sizeof(buffer), "%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", p_cell->STATES[V],
               p_cell->RATES[V], dt, p_cell->STATES[cai] * cml::math::MILLI_TO_NANO, p_cell->ALGEBRAIC[INa] * cml::math::MICRO_TO_NANO,
               p_cell->ALGEBRAIC[INaL] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[ICaL] * cml::math::MICRO_TO_NANO,
               p_cell->ALGEBRAIC[Ito] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[IKr] * cml::math::MICRO_TO_NANO,
               p_cell->ALGEBRAIC[IKs] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[IK1] * cml::math::MICRO_TO_NANO, inet, inet_apd);
      fprintf(fp_time_series, "%.0lf,%s", floor(tcurr), buffer);
      p_features.vm_data.insert(std::pair<double, double>(tcurr, p_cell->STATES[V]));
      p_features.cai_data.insert(std::pair<double, double>(tcurr, p_cell->STATES[cai]));
      next_print_time += dtw;
    }

    icount++;
  }  // end computation loop

  // Assign the remaining features
  // and write it into file.
  snprintf(buffer, sizeof(buffer), "%s/%s/%.2lf/%s_%.2lf_features_core%d_%s.csv", cml::commons::RESULT_FOLDER, p_param->drug_name, conc, p_param->drug_name, conc,
           MPI_Profile::rank, p_param->user_name);
  collect_features(p_features, p_param, p_cell, conc, inet_auc, inet_apd_auc, inal_auc, ical_auc, inal_auc_control, ical_auc_control, buffer,
                   sample_id, group_id);

  if (solver_type == 1) {
    clean_cvode(p_cvode);
    delete p_cvode;
  }
  fclose(fp_time_series);
}

void get_vm_features_postprocessing(Cellmodel *p_cell, Cipa_Features &p_features, const double tcurr) {
  if (tcurr > (p_cell->CONSTANTS[stim_start] + 10) && tcurr < (p_cell->CONSTANTS[stim_start] + 30)) {
    if (p_cell->STATES[V] > p_features.vm_peak) {
      p_features.vm_peak = p_cell->STATES[V];
      if (p_features.vm_peak > 0) {
        p_features.vm_amp30 = p_features.vm_peak - (0.3 * (p_features.vm_peak - p_features.vm_valley));
        p_features.vm_amp50 = p_features.vm_peak - (0.5 * (p_features.vm_peak - p_features.vm_valley));
        p_features.vm_amp90 = p_features.vm_peak - (0.9 * (p_features.vm_peak - p_features.vm_valley));
        p_features.time_vm_peak = tcurr;
      }
    }
  }

  if (p_cell->RATES[V] > p_features.dvmdt_peak) p_features.dvmdt_peak = p_cell->RATES[V];

  if (tcurr > (p_cell->CONSTANTS[stim_start] + 30)) {
    if (p_cell->RATES[V] > p_features.dvmdt_repol_max && p_features.vm_amp30 >= p_cell->STATES[V] && p_cell->STATES[V] >= p_features.vm_amp90) {
      p_features.dvmdt_repol_max = p_cell->RATES[V];
    }
  }
}

void get_ca_features_postprocessing(Cellmodel *p_cell, Cipa_Features &p_features, const double tcurr) {
  if (p_cell->STATES[cai] * cml::math::MILLI_TO_NANO > p_features.ca_peak) {
    p_features.ca_peak = p_cell->STATES[cai] * cml::math::MILLI_TO_NANO;
    p_features.ca_amp30 = p_features.ca_peak - (0.3 * (p_features.ca_peak - p_features.ca_valley));
    p_features.ca_amp50 = p_features.ca_peak - (0.5 * (p_features.ca_peak - p_features.ca_valley));
    p_features.ca_amp90 = p_features.ca_peak - (0.9 * (p_features.ca_peak - p_features.ca_valley));
    p_features.time_ca_peak = tcurr;
  }
}

void get_duration_postprocessing(Cipa_Features &p_features) {
  double time_duration50_prev;
  double time_duration90_prev;

  std::multimap<double, double>::iterator itrmap;
  for (itrmap = p_features.vm_data.begin(); itrmap != p_features.vm_data.end(); itrmap++) {
    if (itrmap->first < p_features.time_vm_peak) {
      if (itrmap->second < p_features.vm_amp50) time_duration50_prev = itrmap->first;
      if (itrmap->second < p_features.vm_amp90) time_duration90_prev = itrmap->first;
    } else {
      if (itrmap->second > p_features.vm_amp50) p_features.apd50 = itrmap->first - time_duration50_prev;
      if (itrmap->second > p_features.vm_amp90) p_features.apd90 = itrmap->first - time_duration90_prev;
    }
  }
  for (itrmap = p_features.cai_data.begin(); itrmap != p_features.cai_data.end(); itrmap++) {
    if (itrmap->first < p_features.time_ca_peak) {
      if (itrmap->second * cml::math::MILLI_TO_NANO < p_features.ca_amp50) time_duration50_prev = itrmap->first;
      if (itrmap->second * cml::math::MILLI_TO_NANO < p_features.ca_amp90) time_duration90_prev = itrmap->first;
    } else {
      if (itrmap->second * cml::math::MILLI_TO_NANO > p_features.ca_amp50) p_features.cad50 = itrmap->first - time_duration50_prev;
      if (itrmap->second * cml::math::MILLI_TO_NANO > p_features.ca_amp90) p_features.cad90 = itrmap->first - time_duration90_prev;
    }
  }
}

void collect_features(Cipa_Features &p_features, const Parameter *p_param, Cellmodel *p_cell, double conc, double inet_auc, double inet_apd_auc,
                      double inal_auc, double ical_auc, double inal_auc_control, double ical_auc_control, char *features_file_name, short sample_id,
                      short group_id) {
  bool file_exists_and_not_empty = false;
  {
    std::ifstream infile(features_file_name);
    file_exists_and_not_empty = infile.good() && infile.peek() != std::ifstream::traits_type::eof();
  }
  char buffer[200];
  FILE *fp_features;

  p_features.qnet = inet_auc / 1000.;
  p_features.qnet_apd = inet_apd_auc / 1000.;
  p_features.vm_valley = p_cell->STATES[V];
  p_features.ca_valley = p_cell->STATES[cai] * cml::math::MILLI_TO_NANO;
  get_duration_postprocessing(p_features);
  p_features.apd_triangulation = p_features.apd90 - p_features.apd50;
  p_features.cad_triangulation = p_features.cad90 - p_features.cad50;
  p_features.inal_auc = inal_auc;
  p_features.ical_auc = ical_auc;

  if ((int)(floor(conc)) == 0) {
    inal_auc_control = inal_auc;
    ical_auc_control = ical_auc;
  }
  p_features.qinward = ((p_features.inal_auc / inal_auc_control) + (p_features.ical_auc / ical_auc_control)) * 0.5;

  fp_features = fopen(features_file_name, "a");
  if (!file_exists_and_not_empty) {
    fprintf(fp_features, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", "sample", "qnet", "qnet_apd", "qInward", "inal_auc", "ical_auc",
            "apd90", "apd50", "apd_tri", "vm_peak", "vm_valley", "dvmdt_peak", "dvmdt_max_repol", "cad90", "cad50", "cad_tri", "ca_peak",
            "ca_valley");
  }
  fprintf(fp_features, "%hd,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", sample_id,
          p_features.qnet, p_features.qnet_apd, p_features.qinward, p_features.inal_auc, p_features.ical_auc, p_features.apd90, p_features.apd50,
          p_features.apd_triangulation, p_features.vm_peak, p_features.vm_valley, p_features.dvmdt_peak, p_features.dvmdt_repol_max, p_features.cad90,
          p_features.cad50, p_features.cad_triangulation, p_features.ca_peak, p_features.ca_valley);

  fclose(fp_features);
}
