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
  snprintf(latex_filename,sizeof(latex_filename),"report_drug_%s_%s_%s.tex",
           p_param->drug_name, p_param->cell_model, p_param->user_name);
  FILE *fp_latex = fopen(latex_filename,"w");
  short font_size = 11;
  float image_scale = 0.8;
  float image_horizontal_padding = 0.;
  const char *DOCUMENT_TITLE = "CardioSim Single Cell Drug Simulation Result";
  const char *USER_NAME = p_param->user_name;

  // header part of LaTEX file
  fprintf(fp_latex,"\\documentclass[%dpt]{article}\n",font_size);
  fprintf(fp_latex,"\\usepackage{graphicx}\n");
  fprintf(fp_latex,"\\usepackage{float}\n");
  fprintf(fp_latex,"\\title{%s}\n",DOCUMENT_TITLE);
  fprintf(fp_latex,"\\author{%s}\n",USER_NAME);
  
  // start of LaTEX document
  fprintf(fp_latex,"\n\\begin{document}\n");
  fprintf(fp_latex,"\\maketitle\n");
  fprintf(fp_latex,"\\section*{Simulation Parameter}\n");
  fprintf(fp_latex,"Cellmodel: %s\n\\\\Drug Name: %s\n\\\\Cmax: %s\n\\\\Ion Channel: 4 (INa, INaL, ICaL, IKr)\n", 
                               (cml::commons::MAP_CELL_NAME.at(p_param->cell_model)).c_str(), p_param->drug_name,p_param->drug_concentrations);

  fprintf(fp_latex,"\\section*{Simulation Protocol}\n");
  fprintf(fp_latex,"Each drug was simulated by inducing %d beats at a %.0lf ms cycle.\n", p_param->number_pacing, p_param->cycle_length);
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
  fprintf(fp_latex,"\\\\(6)  CaD90: Calcium Duration at 90\\%% recovery. Time it takes for the concentration of calcium ions within a heart cell to recover to 90\\%% of its peak level.\n");
  fprintf(fp_latex,"\\\\(7)  CaD50: Calcium Duration at 50\\%% recovery. Time it takes for the concentration of calcium ions within a heart cell to recover to 50\\%% of its peak level\n");
  fprintf(fp_latex,"\\\\(8)  CaDtri: Triangulation of the CaD. Difference between CaD90 and CaD50.\n");

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
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/Vm_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Action potential result in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_vm}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/dVm_dt_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Rates of action potential result in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_dvmdt}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/Cai_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Intracellular calcium concentration result in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_cai}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/INaL_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{L-type sodium current result in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_inal}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/INa_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Sodium current result in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_ina}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/ICaL_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{L-type calcium current result in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_ical}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/IKs_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Slow potassium current result in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_iks}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/IK1_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Inward rectifier current result in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_ik1}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/Ito_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Transient outward current result in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_ito}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/IKr_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{hERG current result in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_ikr}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/%s/Inet_plot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Cummulative current result (iNet) in various concentrations.}\n");
  fprintf(fp_latex,"\t\\label{fig:time_series_inet}\n");
  fprintf(fp_latex,"\\end{figure}\n");


  fprintf(fp_latex,"\\subsection*{Feature Distribution}\n");
  fprintf(fp_latex,"The biomarker for cardiotoxicity assessment proposed by the US FDA Leading Study Group is qNet, and the following figures show the distribution of qNet and other features.\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/%s/qnet_boxplot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Distribution of qNet at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_qnet}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/%s/apd90_boxplot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Distribution of APD90 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_apd90}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/%s/apd50_boxplot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Distribution of APD50 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_apd50}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/%s/cad90_boxplot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Distribution of CaD90 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_cad90}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/%s/cad50_boxplot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Distribution of CaD50 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_cad50}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/%s/vm_peak_boxplot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Distribution of peak action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_vm_peak}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/%s/vm_valley_boxplot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Distribution of lowest action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_vm_valley}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/%s/ca_peak_boxplot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Distribution of peak calcium concentration at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_ca_peak}\n");
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/%s/ca_valley_boxplot.png}\n",
                    image_scale,p_param->drug_name);
  fprintf(fp_latex,"\t\\caption{Distribution of lowest calcium concentrations  at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_ca_valley}\n");
  fprintf(fp_latex,"\\end{figure}\n");


  // end of LaTEX document
  fprintf(fp_latex,"\n\\end{document}\n");
  
  fclose(fp_latex);
/*  
  char command[300];
  snprintf(command, sizeof(command), "pdflatex %s > pdflatex_output.log 2>&1", latex_filename);
  int result = std::system(command);
  if(result == 0){
    printf("PDF generated succesfully!!\n");
  }
  else{
    printf("Failed to generate PDF!!\n");
  }
*/
  return 0;
}
