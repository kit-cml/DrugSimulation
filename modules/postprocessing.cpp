#include "postprocessing.hpp"

#if defined CIPAORDV1_0
#include <types/cellmodels/ohara_rudy_cipa_v1_2017.hpp>
#elif defined TOR_ORD
#include <types/cellmodels/Tomek_model.hpp>
#elif defined TOR_ORD_DYNCL
#include <types/cellmodels/Tomek_dynCl.hpp>
#elif defined GRANDI
#include <types/cellmodels/grandi_2011_atrial_with_meta.hpp>
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

int postprocessing(double conc, double inal_auc_control, double ical_auc_control, const Drug_Row &hill, const Drug_Row &herg,
                    const Parameter *p_param, Cipa_Features &p_features, short sample_id, short group_id, const Cvar_Row *cvar) {
  bool is_ead = false;
  // simulation parameters.
  // Assigned to the constant variables for
  // the sake of simplicity.
  const double cycle_length = p_param->cycle_length;
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
#elif defined GRANDI
  p_cell = new grandi_2011_atrial_with_meta();
  if (is_cvar && cvar) {
    p_cell->initConsts(conc, hill.data, cvar->data);
  } else {
    p_cell->initConsts(conc, hill.data);
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
  FILE *fp_time_series;

  // replace the initial condition
  // with the last state value from
  // the in-silico simulation
  mpi_printf(cml::commons::MASTER_NODE, "STATES before:\n");
  for (short idx = 0; idx < 20; idx++) {
    mpi_printf(cml::commons::MASTER_NODE, "%lf ", p_cell->STATES[idx]);
  }
  mpi_printf(cml::commons::MASTER_NODE, "\nSTATES from initial_values vector:\n");
  for (short idx = 0; idx < 20; idx++) {
    mpi_printf(cml::commons::MASTER_NODE, "%lf ", p_features.initial_values[idx]);
  }
  mpi_printf(cml::commons::MASTER_NODE, "\nUsing initial values from the in-silico simulation.\n");

  if( p_features.initial_values.size() != p_cell->states_size ){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr, "Data mismatch between cell model and initial values file!! Init_Size: %d States_Size_: %d\n",buffer, p_features.initial_values.size(), p_cell->states_size);
    for (int idx = 0; idx < p_features.initial_values.size() ; idx++) {
      mpi_printf(0, "%lf ", p_features.initial_values[idx]);
    }
    mpi_printf(0, "\n");
    for (int idx = 0; idx < p_cell->states_size ; idx++) {
      mpi_printf(0, "%lf ", p_cell->states_size);
    }
    return 1;    
  }
  else {
    copy(p_features.initial_values.begin(), p_features.initial_values.end(), p_cell->STATES);
  }

  for( int idx = 0; idx < p_cell->states_size; idx++ ){
    p_cell->STATES[idx] = p_features.initial_values[idx];
  }
  //set_initial_condition_postprocessing(p_cell, buffer);

  mpi_printf(cml::commons::MASTER_NODE, "STATES after:\n");
  for (short idx = 0; idx < 10; idx++) {
    mpi_printf(cml::commons::MASTER_NODE, "%lf ", p_cell->STATES[idx]);
  }
  mpi_printf(cml::commons::MASTER_NODE, "\n");

  snprintf(buffer, sizeof(buffer), "%s/%.2lf/%s_%.2lf_time_series_smp%d.csv", 
          cml::commons::RESULT_FOLDER, conc, drug_name, conc, sample_id);
  fp_time_series = fopen(buffer, "w");
  if(fp_time_series == NULL){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr, "Cannot create file %s. Make sure the directory is existed!!!\n",buffer);
    return 1;
  }
#if defined GRANDI
  fprintf(fp_time_series, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", 
          "Time(ms)", "Vm(mV)", "dVm/dt(mV/ms)",
          "Cai(nM)", "INa(A/F)", "INaL(mA/F)", "ICaL(A/F)",
          "Ito(A/F)", "IKr(mA/F)", "IKs(mA/F)", "IK1(A/F)", "Inet(A/F)","Inet_AUC(C/F)");
#else
  fprintf(fp_time_series, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", 
          "Time(ms)", "Vm(mV)", "dVm/dt(mV/ms)",
          "Cai(nM)", "INa(nA/uF)", "INaL(nA/uF)", "ICaL(nA/uF)",
          "Ito(nA/uF)", "IKr(nA/uF)", "IKs(nA/uF)", "IK1(nA/uF)", "Inet(uA/uF)","Inet_AUC(uC/uF)");
#endif

#if defined GRANDI
  mpi_printf(0,"gKr after drug in postprocessing: %.4lf\n", p_cell->CONSTANTS[gKr]);
#endif


  // time variables.
  // some of these values are taken from the supplementary materials of ORd2011
  // page 19.
  double time_step = time_step_min;
  double tcurr = 0.;
  double time_point = 25.;
  double tmax = cycle_length;
  short pace_count = 0;
  double next_print_time = 0.0;
  double start_time = 0.0;
  double next_output_time = start_time;
  double tprint = 0.0;

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
    mpi_fprintf(0, stderr, "Solver type %s is undefined! Please choose the available solver types from the manual!\n", solver_type);
    return 1;
  }


  double inal_auc = 0.;
  double ical_auc = 0.;
  double inet = 0.;
  double inet_auc = 0.;
  double inet_apd = 0.;
  double inet_apd_auc = 0.;
  p_features.init(p_cell->STATES[V], p_cell->STATES[cai] * cml::math::MILLI_TO_NANO);
  p_features.vm_amp90 = -77.;
  
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

    // calculating the AUC
    inet = p_cell->ALGEBRAIC[INaL] + p_cell->ALGEBRAIC[ICaL] + p_cell->ALGEBRAIC[Ito] + p_cell->ALGEBRAIC[IKr] + p_cell->ALGEBRAIC[IKs] +
           p_cell->ALGEBRAIC[IK1];
    inet_auc += inet * time_step;
    inal_auc += p_cell->ALGEBRAIC[INaL] * time_step;
    ical_auc += p_cell->ALGEBRAIC[ICaL] * time_step;
    if (p_cell->STATES[V] > p_features.vm_amp90) {
      inet_apd = p_cell->ALGEBRAIC[INaL] + p_cell->ALGEBRAIC[ICaL] + p_cell->ALGEBRAIC[Ito] + p_cell->ALGEBRAIC[IKr] + p_cell->ALGEBRAIC[IKs] +
                 p_cell->ALGEBRAIC[IK1];
      inet_apd_auc += inet_apd * time_step;
    }

    // finding other vm and cai features
    get_vm_features_postprocessing(p_cell, p_features, tcurr);
    get_ca_features_postprocessing(p_cell, p_features, tcurr);
    p_features.vm_data.insert(std::pair<double, double>(tprint, p_cell->STATES[V]));
    p_features.cai_data.insert(std::pair<double, double>(tprint, p_cell->STATES[cai]));

    if( tcurr >= next_output_time - cml::math::EPSILON ){
      tprint = next_output_time-start_time;
#if defined GRANDI
      snprintf(buffer, sizeof(buffer), 
               "%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", 
               p_cell->STATES[V],p_cell->RATES[V], p_cell->STATES[cai] * cml::math::MILLI_TO_NANO, 
               p_cell->ALGEBRAIC[INa], p_cell->ALGEBRAIC[INaL] * cml::math::BASE_TO_MILLI , p_cell->ALGEBRAIC[ICaL],
               p_cell->ALGEBRAIC[Ito], p_cell->ALGEBRAIC[IKr] * cml::math::BASE_TO_MILLI,
               p_cell->ALGEBRAIC[IKs] * cml::math::BASE_TO_MILLI, p_cell->ALGEBRAIC[IK1], inet, inet_auc);
#else
      snprintf(buffer, sizeof(buffer), 
               "%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", p_cell->STATES[V],
               p_cell->RATES[V], p_cell->STATES[cai] * cml::math::MILLI_TO_NANO, p_cell->ALGEBRAIC[INa],
               p_cell->ALGEBRAIC[INaL] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[ICaL] * cml::math::MICRO_TO_NANO,
               p_cell->ALGEBRAIC[Ito] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[IKr] * cml::math::MICRO_TO_NANO,
               p_cell->ALGEBRAIC[IKs] * cml::math::MICRO_TO_NANO, p_cell->ALGEBRAIC[IK1] * cml::math::MICRO_TO_NANO, inet, inet_auc);
#endif
      fprintf(fp_time_series, "%.4lf,%s", tprint, buffer);
      next_output_time += writing_step;
    }

  }  // end computation loop

  // Assign the remaining features
  // and write it into file.
  snprintf(buffer, sizeof(buffer), "%s/%.2lf/%s_%.2lf_features_core%d.csv", 
          cml::commons::RESULT_FOLDER, conc, drug_name, conc,
          MPI_Profile::rank);
  collect_features(p_features, p_param, p_cell, conc, inet_auc, inet_apd_auc, inal_auc, ical_auc, inal_auc_control, ical_auc_control, buffer,
                   sample_id, group_id);

  if (strncasecmp(solver_type, "CVode", sizeof(solver_type)) == 0) {
    clean_cvode(p_cvode);
    delete p_cvode;
  }
  fclose(fp_time_series);

  return 0;
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
                      double inal_auc, double ical_auc, double inal_auc_control, double ical_auc_control, char *features_file_name, 
                      short sample_id, short group_id) 
{

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

  // assuming control concentration value is arbitrary small.
  if (fabs(conc) < cml::math::EPSILON) {
    inal_auc_control = inal_auc;
    ical_auc_control = ical_auc;
  }
  p_features.qinward = ((p_features.inal_auc / inal_auc_control) + (p_features.ical_auc / ical_auc_control)) * 0.5;

  fp_features = fopen(features_file_name, "a");
  if( fp_features == NULL ){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr, "Cannot create file %s. Make sure the directory is existed!!!\n", features_file_name);
    return;
  }
  if (!file_exists_and_not_empty) {
    fprintf(fp_features, "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", "sample", "qNet", "qNetAPD", "qInward", "INaLAUC", "ICaLAUC",
            "APD90", "APD50", "APDtri", "Vmmax", "Vmrest", "dVmdtmax", "dVmdtmaxrepol", "CaTD90", "CaTD50", "CaTDtri", "Camax",
            "Carest");
  }
  fprintf(fp_features, "%hd,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf,%.4lf\n", sample_id,
          p_features.qnet, p_features.qnet_apd, p_features.qinward, p_features.inal_auc, p_features.ical_auc, p_features.apd90, p_features.apd50,
          p_features.apd_triangulation, p_features.vm_peak, p_features.vm_valley, p_features.dvmdt_peak, p_features.dvmdt_repol_max, p_features.cad90,
          p_features.cad50, p_features.cad_triangulation, p_features.ca_peak, p_features.ca_valley);

  fclose(fp_features);
}
