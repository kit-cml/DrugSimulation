#include "report_drug.hpp"

#include <types/cml_consts.hpp>
#include <types/cipa_features.hpp>
#include <types/drug_block_input.hpp>
#include <types/mpi_profile.hpp>
#include <functions/inputoutput.hpp>

#include <cstring>
#include <cstdio>
#include <ctime>
#include <map>
#include <string>
#include <vector>

using std::multimap;
using std::string;
using std::vector;

int generate_report_drug(const Parameter *p_param)
{
  time_t time_data = time(NULL);
  struct tm *current_time = localtime(&time_data);

  char latex_filename[255];
  //snprintf(latex_filename,sizeof(latex_filename),"report_drug_%s_%s_%04d%02d%02d%02d%02d%02d.tex",
  //         p_param->drug_name, p_param->user_name, current_time->tm_year + 1900, current_time->tm_mon + 1, current_time->tm_mday,
  //         current_time->tm_hour, current_time->tm_min, current_time->tm_sec);
  snprintf(latex_filename,sizeof(latex_filename),"./%s/report_drug_%s_%s.tex", cml::commons::RESULT_FOLDER, p_param->drug_name, p_param->cell_model);
  FILE *fp_latex = fopen(latex_filename,"w");
  if( fp_latex == NULL ){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr, "Cannot create file %s. Make sure the directory is existed!!!\n", latex_filename);
    return 1;
  }
  short font_size = 11;
  float image_scale = 0.8;
  float image_horizontal_padding = 0.;
  const char *DOCUMENT_TITLE = "CardioSim Single Cell Drug Simulation Result";
  const char *USER_NAME = p_param->user_name;
  const char *CELL_MODEL_NAME = (cml::commons::MAP_CELL_NAME.at(p_param->cell_model)).c_str();
  const char *DRUG_NAME = p_param->drug_name;
  const char *DRUG_CONCENTRATIONS = p_param->drug_concentrations;
  const int DRUG_CONCENTRATIONS_SIZE = get_concentrations_size(p_param->drug_concentrations);

  // header part of LaTEX file
  fprintf(fp_latex,"\\documentclass[%dpt]{article}\n",font_size);
  fprintf(fp_latex,"\\usepackage[a4paper, top=2.5cm, bottom=2.5cm, left=3cm, right=3cm]{geometry}\n");
  fprintf(fp_latex,"\\usepackage{graphicx}\n");
  fprintf(fp_latex,"\\usepackage{float}\n");
  fprintf(fp_latex,"\\title{%s}\n",DOCUMENT_TITLE);
  fprintf(fp_latex,"\\author{%s}\n",USER_NAME);
  
  // start of LaTEX document
  fprintf(fp_latex,"\n\\begin{document}\n");
  fprintf(fp_latex,"\\maketitle\n");
  fprintf(fp_latex,"\\section*{Simulation Parameter}\n");
  fprintf(fp_latex,"Cellmodel: \\textbf{%s}\n\\\\Drug Name: \\textbf{%s}\n\\\\Cmax: \\textbf{%s}\n\\\\Ion Channel: 4 (INa, INaL, ICaL, IKr)\n", 
                               CELL_MODEL_NAME, DRUG_NAME, DRUG_CONCENTRATIONS);

  fprintf(fp_latex,"\\section*{Simulation Protocol}\n");
  fprintf(fp_latex,"Each drug was simulated by inducing \\textbf{%d} beats at a \\textbf{%.0lf} ms cycle.\n", p_param->number_pacing, p_param->cycle_length);
  fprintf(fp_latex,"\\\\Biomarker: qNet, qInward, APD90, APD50, APDtri, CaD90, CaD50, CaDtri.\n");
  fprintf(fp_latex,"\\\\(1) qNet: Net charge carried by six major currents over a simulated beat (IKr, ICaL, INaL, Ito, IKs, IK1).\n");
  fprintf(fp_latex,"\\begin{equation}\n");
  fprintf(fp_latex,"\t qNet = \\int_{0}^{CL} (I_{NaL} + I_{CaL} + I_{Kr} + I_{K1} + I_{Ks} + I_{to}) \\,dt\n");
  fprintf(fp_latex,"\t \\label{qnet}\n");
  fprintf(fp_latex,"\\end{equation}\n");
  fprintf(fp_latex,"\\\\(2) qInward: Inward currents carried by cations (such as sodium or calcium) into the cell during the depolarization phase of the cardiac action potential.\n");
  fprintf(fp_latex,"\\\\(3) APD90: Action Potential Duration at 90\\%% repolarization. time it takes for the action potential of a cardiac cell to repolarize to 90\\%% of its maximum\n");
  fprintf(fp_latex,"\\\\(4) APD50: APD at 50\\%% repolarization. time it takes for the action potential of a cardiac cell to repolarize to 50\\%% of its maximum\n");
  fprintf(fp_latex,"\\\\(5)  APDtri: Triangulation of the APD. Difference between APD90 and APD50. \n");
  fprintf(fp_latex,"\\\\(6)  CaTD90: Calcium Transient Duration at 90\\%% recovery. Time it takes for the concentration of calcium ions within a heart cell to recover to 90\\%% of its peak level.\n");
  fprintf(fp_latex,"\\\\(7)  CaTD50: Calcium Transient Duration at 50\\%% recovery. Time it takes for the concentration of calcium ions within a heart cell to recover to 50\\%% of its peak level\n");
  fprintf(fp_latex,"\\\\(8)  CaTDtri: Triangulation of the CaTD. Difference between CaD90 and CaD50.\n");

  fprintf(fp_latex,"\\section*{Method}\n");
  fprintf(fp_latex,"Calculate biomarkers by simulating one cycle after assigning the state values of the ventricular cell membrane protein gates derived through in silico simulation as the initial conditions of the model.\n");


  fprintf(fp_latex,"\\section*{Result}\n");
  fprintf(fp_latex,"\\subsection*{Time-series Result}\n");
  fprintf(fp_latex,"The following results will show the time-series result of the action potential, rates of action potentials, intracellular calcium, and ionic currents involved in the drugs.\n");


  // generate time series for control and other concentrations.
  // Select 5 of random samples to be displayed.
  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/Vm_plot.png}\n",image_scale);
  fprintf(fp_latex,"\t\\caption{Action potential profiles of the \\textbf{%s} model at \\textbf{%d} \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control condition, while colored solid lines correspond to drug-related simulations.}\n", CELL_MODEL_NAME, DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS );
  fprintf(fp_latex,"\t\\label{fig:time_series_vm}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/dVm_dt_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Time derivative membrane voltage (dV/dt) at various \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_dvmdt}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/Cai_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Intracellular calcium concentrations (\\([Ca^{2+}]_i\\)) traces simulated at \\textbf{%d} \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line indicates the control condition, and the solid lines represent drug-related simulations }\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_cai}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/INaL_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Late sodium current (INaL) traces at different \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_inal}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/INa_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Fast sodium current (INa) traces at \\textbf{%d} \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_ina}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/ICaL_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{L-type calcium current (ICaL) traces at \\textbf{%d} \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_ical}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/IKs_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Slow delayed rectifier potassium current (IKs) at different \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_iks}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/IK1_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Inward rectifier potassium current (IK1) at various \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_ik1}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/Ito_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Transient outward current (Ito) at various \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_ito}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/IKr_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{hERG current (IKr) at \\textbf{%d} \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_ikr}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/Inet_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Net current between inward and outward channels (iNet) at \\textbf{%d} \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_inet}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/Inet_AUC_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{AUC of iNet in time-series data (qNet) at various \\textbf{%s} concentrations (\\textbf{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_inet_auc}\n");
  fprintf(fp_latex,"\\end{figure}\n");


  fprintf(fp_latex,"\\subsection*{Feature Distribution}\n");
  fprintf(fp_latex,"The biomarker for cardiotoxicity assessment proposed by the US FDA Leading Study Group is qNet, and the following figures show the distribution of qNet and other features.\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/qNet_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of qNet at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_qnet}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/qInward_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of qInward at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_qinward}\n");
  fprintf(fp_latex,"\\end{figure}\n");


  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/APD90_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of APD90 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_apd90}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/APD50_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of APD50 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_apd50}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/APD_tri_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of APD triangulation at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_apdtri}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/Vm_max_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_vm_max}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/Vm_rest_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of minimum action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_vm_rest}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/dVmdt_max_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum rates action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_dvmdt_max}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/dVmdt_max_repol_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum rates action potential during repolarization at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_dvmdt_max_repol}\n");
  fprintf(fp_latex,"\\end{figure}\n");


  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTD90_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of CaTD90 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_catd90}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTD50_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of CaTD50 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_catd50}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTD_tri_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of CaTD triangulation at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_catdtri}\n");
  fprintf(fp_latex,"\\end{figure}\n");


  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/Ca_max_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum calcium concentration at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_ca_max}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/Ca_rest_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of minimum calcium concentrations  at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_ca_valley}\n");
  fprintf(fp_latex,"\\end{figure}\n");


  // end of LaTEX document
  fprintf(fp_latex,"\n\\end{document}\n");
  
  fclose(fp_latex);
  return 0;
}

int get_concentrations_size(const char* concs_str)
{
  int count = 1;
  if (concs_str == 0) count = 0;
  else{
    for (const char* pdx = concs_str; *pdx; ++pdx)
        if (*pdx == ',') count++;
  }

  return count;
}
