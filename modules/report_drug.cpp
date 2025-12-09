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
  fprintf(fp_latex,"\\usepackage{bm}\n");
  fprintf(fp_latex,"\\usepackage{booktabs}\n");
  fprintf(fp_latex,"\\usepackage{graphicx}\n");
  fprintf(fp_latex,"\\usepackage{float}\n");
  fprintf(fp_latex,"\\usepackage{pgfplotstable}\n");
  fprintf(fp_latex,"\\usepackage{xcolor}\n");
  fprintf(fp_latex,"\\title{%s}\n",DOCUMENT_TITLE);
  fprintf(fp_latex,"\\author{%s}\n",USER_NAME);
  
  // start of LaTEX document
  fprintf(fp_latex,"\n\\begin{document}\n");
  fprintf(fp_latex,"\\maketitle\n");
#if defined teatrseatsetaet
  fprintf(fp_latex,"\\section*{Simulation Parameter}\n");
  fprintf(fp_latex,"Cellmodel: \\textcolor{blue}{%s}\n\\\\Drug Name: \\textcolor{blue}{%s}\n\\\\Concentrations (nM): \\textcolor{blue}{%s}\n\\\\Ion Channels: 4 (INa, INaL, ICaL, IKr)\n", 
                               CELL_MODEL_NAME, DRUG_NAME, DRUG_CONCENTRATIONS);

  fprintf(fp_latex,"\\section*{Method}\n");
#endif
  fprintf(fp_latex,"\\subsection*{Cell Model and Simulation Environment}\n");
  fprintf(fp_latex,"Simulations were performed using the human atrial myocyte model developed by \\textcolor{blue}{%s}, which incorporates detailed representations of major ionic currents (I$_{Na}$, I$_{NaL}$, I$_{CaL}$, I$_{to}$, I$_{Kr}$, I$_{Ks}$, I$_{K1}$), ion transporters such as the Na$^{+}$/Ca$^{2+}$ exchanger (NCX) and Na$^{+}$/K$^{+}$ ATPase (NKA), the SERCA pump, and intracellular Ca$^{2+}$ handling processes involving the junctional and network sarcoplasmic reticulum (JSR and NSR).\n\n", CELL_MODEL_NAME);
  fprintf(fp_latex,"All differential equations were solved using a stiff ODE solver (CVODE, Backward Differentiation Formula). The basic integration time step was set to \\textcolor{blue}{%.2lf} ms to ensure numerical stability and precision.\n", TIME_STEP);

  fprintf(fp_latex,"\\subsection*{Drug Modeling and Application}\n");
  fprintf(fp_latex,"\\textcolor{blue}{%s} was used as the test compound. Its inhibitory effects on ionic currents were implemented by modifying the maximal conductance (g$_{X}$) of each affected current using the Hill equation. For each target current (I$_{Na}$, I$_{NaL}$, I$_{CaL}$, I$_{Kr}$), \\textbf{$\\bm{IC_{50}}$ values and Hill coefficients (h) estimated from the bootstrap module} were applied:.\n", DRUG_NAME);
  fprintf(fp_latex,"\\begin{equation}\n");
  fprintf(fp_latex,"\t Block = \\frac{1}{1 + \\left( \\frac{C}{IC_{50}}\\right)^{h}}\n");
  fprintf(fp_latex,"\t \\label{drug_block}\n");
  fprintf(fp_latex,"\\end{equation}\n");
  fprintf(fp_latex,"\\\\The drug-modified conductance was then computed as:\n");
  fprintf(fp_latex,"\\begin{equation}\n");
  fprintf(fp_latex,"\t g_{X,drug} = g_{X,control} \\cdot  (1-Block)\n");
  fprintf(fp_latex,"\t \\label{conductance_drug_block}\n");
  fprintf(fp_latex,"\\end{equation}\n");
  fprintf(fp_latex,"\\\\The drug concentration was set to the reported maximal plasma concentration (Cmax = \\textcolor{blue}{%s}) and maintained constant throughout the simulation.\n", DRUG_CONCENTRATIONS);

  fprintf(fp_latex,"\\subsection*{Simulation Protocol}\n");
  fprintf(fp_latex,"For each drug condition, the model was paced at a basic cycle length (BCL) of \\textcolor{blue}{%.0lf} ms. To ensure that ionic concentrations (Ca$_{i}$, Na$_{i}$, K$_{i}$), gating variables, and membrane potential fully adapted to the presence of quinidine, the cell was stimulated for \\textbf{\\textcolor{blue}{%hd} consecutive beats}. This long-term pacing allowed the system to reach a new steady-state under drug exposure.\n\n", CYCLE_LENGTH, NUMBER_PACING);
  fprintf(fp_latex,"During the simulation, the model monitored the most recent 250 pacing cycles and identified the beat with the largest membrane potential derivative during the repolarization phase. At the conclusion of the simulation, all state variables associated with this selected beat—including membrane potential (Vm), gating variables, junctional and network sarcoplasmic reticulum (JSR/NSR) Ca$^{2+}$ content, and intracellular ion concentrations—were saved and used as the \\textbf{initial conditions} for an additional single-beat simulation. Electrophysiological biomarkers were then computed from this final beat to eliminate transient effects and ensure that the results reflected the steady-state response to the drug .\n");

  fprintf(fp_latex,"\\subsection*{Biomarker Computation}\n");
  fprintf(fp_latex,"\\begin{enumerate}\n");


  fprintf(fp_latex,"\t\\item \\textbf{qNet (Net Charge Transfer)} \nqNet was calculated as the time integral of six major currents (I$_{Kr}$, I$_{CaL}$, I$_{NaL}$, I$_{to}$, I$_{Ks}$, I$_{K1}$):\n"); 
  fprintf(fp_latex,"\\begin{equation}\n");
  fprintf(fp_latex,"\t qNet = \\int_{0}^{CL} (I_{NaL} + I_{CaL} + I_{Kr} + I_{K1} + I_{Ks} + I_{to}) \\,dt\n");
  fprintf(fp_latex,"\t \\label{qnet}\n");
  fprintf(fp_latex,"\\end{equation}\n");
  fprintf(fp_latex,"\tThis metric reflects the overall balance between inward and outward currents throughout the action potential.\n");

  fprintf(fp_latex,"\t\\item \\textbf{qInward (Total Inward Cation Change)} \nqInward was obtained by integrating all inward components of cationic currents during the cardiac cycle:\n"); 
  fprintf(fp_latex,"\\begin{equation}\n");
  fprintf(fp_latex,"\t qInward = \n");
  fprintf(fp_latex,"\t 0.5 \\cdot \\left(\n"); 
  fprintf(fp_latex,"\t\t \\frac{\\int_{0}^{CL} I_{\\mathrm{NaL,drug}}\\, dt}{\\int_{0}^{CL} I_{\\mathrm{NaL,control}}\\, dt}\n\t+\n");
  fprintf(fp_latex,"\t\t \\frac{\\int_{0}^{CL} I_{\\mathrm{CaL,drug}}\\, dt}{\\int_{0}^{CL} I_{\\mathrm{CaL,control}}\\, dt}\n");
  fprintf(fp_latex,"\t \\right)\n");
  fprintf(fp_latex,"\t \\label{qInward}\n");
  fprintf(fp_latex,"\\end{equation}\n");
  fprintf(fp_latex,"\tIt indicates the magnitude of Na$^{+}$ and Ca$^{2+}$ influx and is closely associated with intracellular Ca$^{2+}$ loading.\n");

  fprintf(fp_latex,"\t\\item \\textbf{APD$_{90}$ and APD$_{50}$} \nAPD$_{90}$ and APD$_{50}$ were defined as the time from the upstroke of the action potential to the point at which the membrane potential repolarized to 90\\%% and 50\\%% of the difference between peak and diastolic voltage, respectively.\n"); 

  fprintf(fp_latex,"\t\\item \\textbf{APD$_{Tri}$ (APD Triangulation)} \n"); 
  fprintf(fp_latex,"\\begin{equation}\n");
  fprintf(fp_latex,"\t APD_{tri} = APD_{90} - APD_{50}\n");
  fprintf(fp_latex,"\t \\label{apdtri}\n");
  fprintf(fp_latex,"\\end{equation}\n");
  fprintf(fp_latex,"\tThis index reflects changes in the shape (triangulation) of the action potential and is known to be associated with increased proarrhythmic risk.\n");

 fprintf(fp_latex,"\t\\item \\textbf{CaTD$_{90}$ and CaTD$_{50}$} \n These metrics represent the time from the peak of the Ca$^{2+}$ transient to 90\\%% and 50\\%% recovery toward baseline, indicating the kinetics of Ca$^{2+}$ reuptake and extrusion (SERCA activity, NCX function).\n"); 

  fprintf(fp_latex,"\t\\item \\textbf{CaTDTri (Calcium Transient Triangulation)} \n"); 
  fprintf(fp_latex,"\\begin{equation}\n");
  fprintf(fp_latex,"\t CaTD_{tri} = CaTD_{90} - CaTD_{50}\n");
  fprintf(fp_latex,"\t \\label{catdtri}\n");
  fprintf(fp_latex,"\\end{equation}\n");
  fprintf(fp_latex,"\tCaTDtri evaluates the extent of prolonged intracellular Ca$^{2+}$ elevation, which is related to Ca$^{2+}$ overload and triggered arrhythmogenic events.\n");

  fprintf(fp_latex,"\\end{enumerate}\n");

  fprintf(fp_latex,"\\section*{Results}\n");
  fprintf(fp_latex,"\\subsection*{Time-series}\n");
  fprintf(fp_latex,"Figure 1 displays membrane potential profiles simulated with the \\textcolor{blue}{%s}, elucidating the effects of \\textcolor{blue}{%s} at multiple concentrations (\\textcolor{blue}{%s}) nM. In these traces, the control condition (dashed lines) demonstrates a more rapid repolarization compared to drug-treated samples (solid colored lines). By increasing the \\textcolor{blue}{%s} concentration, there is a clear prolongation of action potential duration and delayed membrane potential recovery, which aligns with quinidine’s established pharmacological impact on cardiac tissues.\n\n", CELL_MODEL_NAME, DRUG_NAME, DRUG_CONCENTRATIONS, DRUG_NAME);
  fprintf(fp_latex,"Building on these findings, Figure 2 illustrates the corresponding time derivatives of the membrane potential ($dV/dt$), thereby facilitating a more detailed assessment of depolarization kinetics under the same experimental scenarios. In both control and \\textcolor{blue}{%s}-exposed samples, pronounced peaks in $dV/dt$ are observed during the action potential upstroke, indicating robust initial depolarization. Notably, \\textcolor{blue}{%s} administration does not considerably affect the maximal upstroke velocity. Thus, the principal influence of quinidine lies in modulating action potential duration and repolarization dynamics rather than the rate of membrane potential rise.\n\n", DRUG_NAME, DRUG_NAME);
  fprintf(fp_latex,"In addition to changes in membrane potential, Figure 3 extends the analysis by presenting simulated traces of intracellular calcium concentration ([$Ca^{2+}_{i}$]) for the same \\textcolor{blue}{%s} concentrations. Similar peak [$Ca^{2+}_{i}$] values are evident following stimulation across all experimental groups, as indicated by the dashed (control) and solid (treated) lines. Importantly, the presence of \\textcolor{blue}{%s} does not result in significant differences in the decay phase or peak calcium levels. As a result, these data suggest that \\textcolor{blue}{%s}, within the context of these parameters, predominantly affects electrical properties of ventricular myocytes while minimally influencing calcium handling.\n\n", DRUG_NAME, DRUG_NAME, DRUG_NAME);
  fprintf(fp_latex,"Furthermore, Figure 4 provides traces of the late sodium current (I$_{NaL}$), highlighting marked suppression in response to \\textcolor{blue}{%s}, particularly at higher concentrations. Similarly, Figure 5 demonstrates that peak fast sodium current (I$_{Na}$) is notably reduced in drug-treated samples, reinforcing \\textcolor{blue}{%s}’s efficacy in sodium channel inhibition and arrhythmia suppression.\n\n", DRUG_NAME, DRUG_NAME);
  fprintf(fp_latex,"Taking these findings together with those for the L-type calcium current, Figure 6 shows that \\textcolor{blue}{%s}’s effects are mostly centered upon sodium channel activity, with only modest reduction observed in peak I$_{CaL}$ responses. This suggests preserved calcium influx during excitation.\n\n", DRUG_NAME);
  fprintf(fp_latex,"Additionally, Figure 7 presents the slow delayed rectifier potassium current (I$_{Ks}$), indicating that \\textcolor{blue}{%s} exposure leads to enhanced peak (I$_{Ks}$) compared to control, a compensatory ion channel response within the repolarization phase.\n\n", DRUG_NAME);
  fprintf(fp_latex,"Next, Figure 8 and Figure 9 provide data on inward rectifier (I$_{K1}$) and transient outward (I$_{to}$) potassium currents, respectively, demonstrating that while \\textcolor{blue}{%s} has a modest depressive effect on inward rectifier current, its impact on the transient outward current is minimal. In Figure 10, prominent inhib\n\n", DRUG_NAME);
  fprintf(fp_latex,"Figure 11 summarizes the cumulative electrophysiological impact of \\textcolor{blue}{%s} by displaying net transmembrane current (I$_{net}$), revealing a reduction in peak amplitudes and overall current flow dynamics in drug-treated samples, which reflects the summed effects of suppressed inward and repolarizing currents across the cardiac cycle\n\n", DRUG_NAME);
  fprintf(fp_latex,"Finally, Figure 12 integrates these individual findings through an examination of the area under the curve (AUC) of net current (q$_{Net}$) over time. In both control and drug-related simulations, q$_{Net}$ traces highlight the cumulative differences in charge movement during depolarization and repolarization. \\textcolor{blue}{%s} exposure at both concentrations results in more negative AUC values, and a greater divergence from the control profile, which quantifies the extended action potential duration and the dampened total ionic activity. These integrated data provide quantitative evidence for \\textcolor{blue}{%s}’s concentration-dependent effects on cardiac electrophysiology, complementing the mechanistic insights from individual channel analyses.\n\n", DRUG_NAME, DRUG_NAME);
  fprintf(fp_latex,"\\subsection*{Features Distribution}\n");
  fprintf(fp_latex,"The biomarker for cardiotoxicity assessment proposed by the US FDA Leading Study Group is qNet, and the following figures show the distribution of qNet and other features.\n");


  fprintf(fp_latex,"\\clearpage\n");

  // generate time series for control and other concentrations.
  // Select 5 of random samples to be displayed.
  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/Vm_plot.png}\n",image_scale);
  fprintf(fp_latex,"\t\\caption{Action potential profiles of the \\textcolor{blue}{%s} model at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, while colored solid lines correspond to drug-related simulations.}\n", CELL_MODEL_NAME, DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS );
  fprintf(fp_latex,"\t\\label{fig:time_series_vm}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/dVm_dt_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Time derivative membrane voltage (dV/dt) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_dvmdt}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/Cai_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Intracellular calcium concentrations (\\([Ca^{2+}]_i\\)) traces simulated at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line indicates the control condition, and the solid lines represent drug-related simulations }\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_cai}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/INaL_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Late sodium current (INaL) traces at different \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_inal}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/INa_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Fast sodium current (INa) traces at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_ina}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/ICaL_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{L-type calcium current (ICaL) traces at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_ical}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/IKs_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Slow delayed rectifier potassium current (IKs) at different \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_iks}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/IK1_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Inward rectifier potassium current (IK1) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_ik1}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/Ito_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Transient outward current (Ito) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_ito}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/IKr_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{hERG current (IKr) at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_ikr}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/Inet_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Net current between inward and outward channels (iNet) at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control condition, and the colored solid lines indicate drug-treated simulations.}\n", DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_inet}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/time_series/Inet_AUC_plot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{AUC of iNet in time-series data (qNet) at various \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM). The dashed line represents the control, and the solid lines denote drug-related simulations.}\n", DRUG_NAME, DRUG_CONCENTRATIONS);
  fprintf(fp_latex,"\t\\label{fig:time_series_inet_auc}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/qNet_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of qNet at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_qnet}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/qInward_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of qInward at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_qinward}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");


  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/APD90_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of APD90 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_apd90}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/APD50_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of APD50 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_apd50}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/APDtri_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of APD triangulation at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_apdtri}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/Vmmax_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_vm_max}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/Vmrest_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of minimum action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_vm_rest}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/dVmdtmax_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum rates action potential at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_dvmdt_max}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/dVmdtmaxrepol_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum rates action potential during repolarization at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_dvmdt_max_repol}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTD90_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of CaTD90 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_catd90}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTD50_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of CaTD50 at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_catd50}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/CaTDtri_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of CaTD triangulation at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_catdtri}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");


  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/Camax_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of maximum calcium concentration at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_ca_max}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\hspace{%.2lfcm}\n", image_horizontal_padding);
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/features/Carest_boxplot.png}\n", image_scale);
  fprintf(fp_latex,"\t\\caption{Distribution of minimum calcium concentrations  at different concentration levels, from control, cmax1, cmax2, cmax3, and cmax4.}\n");
  fprintf(fp_latex,"\t\\label{fig:features_ca_valley}\n");
  fprintf(fp_latex,"\\end{figure}\n");
  fprintf(fp_latex,"\\clearpage\n");

  // table always be the last
  // and caption should be above.
  fprintf(fp_latex,"\n\\begin{table}[ht]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\caption{Average feature values of the \\textcolor{blue}{%s} model at \\textcolor{blue}{%d} \\textcolor{blue}{%s} concentrations (\\textcolor{blue}{%s} nM).}\n", CELL_MODEL_NAME, DRUG_CONCENTRATIONS_SIZE, DRUG_NAME, DRUG_CONCENTRATIONS );
  fprintf(fp_latex,"\t\\label{tab:dose_averages}\n");
  fprintf(fp_latex,"\t\\pgfplotstabletypeset[\n");
  fprintf(fp_latex,"\t\tcol sep=comma,\n");
  fprintf(fp_latex,"\t\theader=true,\n");
  fprintf(fp_latex,"\t\tstring type,\n");
  fprintf(fp_latex,"\t\tcolumns/.style={string type},\n");
  fprintf(fp_latex,"\t\tcolumns/*/.style={string type},\n");
  fprintf(fp_latex,"\t\tevery head row/.style={before row=\\toprule, after row=\\midrule},\n");
  fprintf(fp_latex,"\t\tevery last row/.style={after row=\\bottomrule},\n");
  fprintf(fp_latex,"\t ]{./dose_averages_transposed.csv}\n");
  fprintf(fp_latex,"\\end{table}\n\n");
  fprintf(fp_latex,"\\clearpage\n");


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
