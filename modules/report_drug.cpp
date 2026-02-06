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
  snprintf(latex_filename,sizeof(latex_filename),"./%s/report_template.tex", cml::commons::RESULT_FOLDER);
  FILE *fp_latex = fopen(latex_filename,"w");
  if( fp_latex == NULL ){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr, "Cannot create file %s. Make sure the directory is existed!!!\n", latex_filename);
    return 1;
  }
  snprintf(latex_filename,sizeof(latex_filename),"./%s/figures.tex", cml::commons::RESULT_FOLDER);
  FILE *fp_figures = fopen(latex_filename,"w");
  if( fp_figures == NULL ){
    mpi_fprintf(cml::commons::MASTER_NODE, stderr, "Cannot create file %s. Make sure the directory is existed!!!\n", latex_filename);
    return 1;
  }

  short font_size = 11;
  float image_scale = 0.8;
  float image_horizontal_padding = 0.;
  float table_horizontal_padding = -15.;
  const char *DOCUMENT_TITLE = "CardioSim Report";
  const char *USER_NAME = p_param->user_name;
  const char *CELL_MODEL_NAME = (cml::commons::MAP_CELL_NAME.at(p_param->cell_model)).c_str();
  const char *DRUG_NAME = p_param->drug_name;
  const char *DRUG_CONCENTRATIONS = p_param->drug_concentrations;
  const int DRUG_CONCENTRATIONS_SIZE = get_concentrations_size(p_param->drug_concentrations);
  const short NUMBER_PACING = p_param->number_pacing;
  const double CYCLE_LENGTH = p_param->cycle_length;
  const double TIME_STEP = p_param->time_step_min;

  // header part of LaTEX file
  fprintf(fp_latex,"\\documentclass[%dpt]{article}\n",font_size);
  fprintf(fp_latex,"\\usepackage[a4paper, top=2cm, bottom=1.5cm, left=2cm, right=2cm]{geometry}\n");
  fprintf(fp_latex,"\\usepackage{amsmath}\n");
  fprintf(fp_latex,"\\usepackage{bm}\n");
  fprintf(fp_latex,"\\usepackage{booktabs}\n");
  fprintf(fp_latex,"\\usepackage[T1]{fontenc}\n");
  fprintf(fp_latex,"\\usepackage{graphicx}\n");
  fprintf(fp_latex,"\\usepackage{hyperref}\n");
  fprintf(fp_latex,"\\usepackage[utf8]{inputenc}\n");
  fprintf(fp_latex,"\\usepackage{float}\n");
  fprintf(fp_latex,"\\usepackage{pgfplotstable}\n");
  fprintf(fp_latex,"\\usepackage[strings]{underscore}\n");
  fprintf(fp_latex,"\\usepackage{siunitx}\n");
  fprintf(fp_latex,"\\usepackage{xcolor}\n\n");
  fprintf(fp_latex,"\\DeclareUnicodeCharacter{223C}{\\textasciitilde}\n");
  fprintf(fp_latex,"\\DeclareUnicodeCharacter{2212}{-}\n\n");
  fprintf(fp_latex,"\\title{%s}\n",DOCUMENT_TITLE);
  fprintf(fp_latex,"\\author{%s}\n",USER_NAME);

  // start of LaTEX document
  fprintf(fp_latex,"\n\\begin{document}\n");
  fprintf(fp_figures,"\n\\setlength{\\baselineskip}{0.5cm}\n");
  fprintf(fp_figures,"\n\\setlength{\\parskip}{0.5cm}\n\n");

  fprintf(fp_latex,"\\maketitle\n");
  fprintf(fp_latex,"\\section*{Method}\n");

  fprintf(fp_latex,"\\section*{Results}\n");

  fprintf(fp_latex,"\\section*{Discussions}\n");


  fprintf(fp_latex,"\\clearpage\n");

  // generate time series for control and other concentrations.
  // Select 5 of random samples to be displayed.
  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/Vm-plot.png}\n",image_scale);
  fprintf(fp_latex,"\t\\caption{Action potential profiles of the \\textcolor{blue}{%s} model at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, while colored solid lines correspond to drug-related simulations.}\n", CELL_MODEL_NAME, DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS );
  fprintf(fp_latex,"\t\\label{fig:time-series-vm}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/dVmdt-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Time derivative membrane voltage (dV/dt) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-dvmdt}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/Cai-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Intracellular calcium concentrations (\\([Ca^{2+}]_i\\)) traces simulated at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line indicates the control condition, and the solid lines represent drug-related simulations }\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-cai}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/INaL-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Late sodium current (INaL) traces at different \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-inal}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/INa-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Fast sodium current (INa) traces at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-ina}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/ICaL-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{L-type calcium current (ICaL) traces at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-ical}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/IKs-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Slow delayed rectifier potassium current (IKs) at different \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-iks}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/IK1-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Inward rectifier potassium current (IK1) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-ik1}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/Ito-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Transient outward current (Ito) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-ito}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/IKr-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{hERG current (IKr) at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-ikr}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/Inet-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Net current between inward and outward channels (iNet) at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-inet}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time-series/InetAUC-plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{AUC of iNet in time-series data (qNet) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time-series-inet-auc}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/qNet-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of qNet at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-qnet}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/qInward-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of qInward at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-qinward}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");


  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/APD90-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of APD90 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-apd90}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/APD50-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of APD50 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-apd50}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/APDTri-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of APD triangulation at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-apdtri}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/VmMax-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-vmmax}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/VmRest-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of minimum action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-vmrest}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/dVmdtMax-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum rates action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-dvmdtmax}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/dVmdtMaxRepol-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum rates action potential during repolarization at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-dvmdtmaxrepol}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTD90-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of CaTD90 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-catd90}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTD50-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of CaTD50 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-catd50}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTDTri-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of CaTD triangulation at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-catdtri}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");


  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaMax-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum calcium concentration at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-camax}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaRest-boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of minimum calcium concentrations  at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features-carest}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  // table always be the last
  // and caption should be above.
  fprintf(fp_latex,"\n\\begin{table}[ht]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\caption{Median values and Confidence Interval 95\\%% of of averaged samples obtained from the \\textcolor{blue}{%s} model under \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM).}\n", CELL_MODEL_NAME, DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS );
  fprintf(fp_latex,"\t\\label{tab:median-ci}\n");
  fprintf(fp_latex,"\t\\pgfplotstabletypeset[\n");
  fprintf(fp_latex,"\t\tcol sep=comma,\n");
  fprintf(fp_latex,"\t\theader=true,\n");
  fprintf(fp_latex,"\t\tstring type,\n");
  fprintf(fp_latex,"\t\tcolumns/.style={string type},\n");
  fprintf(fp_latex,"\t\tevery head row/.style={before row=\\toprule, after row=\\midrule},\n");
  fprintf(fp_latex,"\t\tevery last row/.style={after row=\\bottomrule},\n");
  fprintf(fp_latex,"\t ]{./%s-median-ci.csv}\n", DRUG_NAME);
  fprintf(fp_latex,"\\end{table}\n\n");
  fprintf(fp_latex,"\\clearpage\n");


  // end of LaTEX document
  fprintf(fp_latex,"\n\\end{document}\n");

  fprintf(fp_figures,"\\clearpage\n");
  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/Vm-plot.png}\n",image_scale);
  fprintf(fp_figures,"\t\\caption{Action potential profiles of the \\textcolor{blue}{%s} model at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, while colored solid lines correspond to drug-related simulations.}\n", CELL_MODEL_NAME, DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS );
  fprintf(fp_figures,"\t\\label{fig:time-series-vm}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/dVmdt-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Time derivative membrane voltage (dV/dt) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-dvmdt}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/Cai-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Intracellular calcium concentrations (\\([Ca^{2+}]_i\\)) traces simulated at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line indicates the control condition, and the solid lines represent drug-related simulations }\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-cai}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/INaL-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Late sodium current (INaL) traces at different \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-inal}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/INa-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Fast sodium current (INa) traces at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-ina}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/ICaL-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{L-type calcium current (ICaL) traces at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-ical}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/IKs-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Slow delayed rectifier potassium current (IKs) at different \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-iks}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/IK1-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Inward rectifier potassium current (IK1) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-ik1}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/Ito-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Transient outward current (Ito) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-ito}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/IKr-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{hERG current (IKr) at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-ikr}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/Inet-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Net current between inward and outward channels (iNet) at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-inet}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/time-series/InetAUC-plot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{AUC of iNet in time-series data (qNet) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_figures,"\t\\label{fig:time-series-inet-auc}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/qNet-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of qNet at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-qnet}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/qInward-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of qInward at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-qinward}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");


  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/APD90-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of APD90 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-apd90}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/APD50-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of APD50 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-apd50}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/APDTri-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of APD triangulation at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-apdtri}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/VmMax-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of maximum action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-vmmax}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/VmRest-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of minimum action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-vmrest}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/dVmdtMax-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of maximum rates action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-dvmdtmax}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/dVmdtMaxRepol-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of maximum rates action potential during repolarization at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-dvmdtmaxrepol}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTD90-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of CaTD90 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-catd90}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTD50-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of CaTD50 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-catd50}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTDTri-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of CaTD triangulation at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-catdtri}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");


  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/CaMax-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of maximum calcium concentration at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-camax}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\\begin{figure}[H]\n");
  fprintf(fp_figures,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\includegraphics[scale=%.1f]{./plots/features/CaRest-boxplot.png}\n", image_scale);
  fprintf(fp_figures,"\t\\caption{Distribution of minimum calcium concentrations  at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_figures,"\t\\label{fig:features-carest}\n");
  fprintf(fp_figures,"\\end{figure}\n");
  fprintf(fp_figures,"\\clearpage\n");

  fprintf(fp_figures,"\n\\begin{table}[ht]\n");
  fprintf(fp_figures,"\t\\centering\n");
  fprintf(fp_figures,"\t\\caption{Median values and Confidence Interval 95\\%% of the \\textcolor{blue}{%s} model at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM).}\n", CELL_MODEL_NAME, DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS );
  fprintf(fp_figures,"\t\\label{tab:median-ci}\n");
  fprintf(fp_figures,"\t\\pgfplotstabletypeset[\n");
  fprintf(fp_figures,"\t\tcol sep=comma,\n");
  fprintf(fp_figures,"\t\theader=true,\n");
  fprintf(fp_figures,"\t\tstring type,\n");
  fprintf(fp_figures,"\t\tcolumns/.style={string type},\n");
  fprintf(fp_figures,"\t\tevery head row/.style={before row=\\toprule, after row=\\midrule},\n");
  fprintf(fp_figures,"\t\tevery last row/.style={after row=\\bottomrule},\n");
  fprintf(fp_figures,"\t ]{./%s-median-ci.csv}\n", DRUG_NAME);
  fprintf(fp_figures,"\\end{table}\n\n");
  fprintf(fp_figures,"\\clearpage\n");

  
  fclose(fp_latex);
  fclose(fp_figures);
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
