#include "report_drug.hpp"

#include <types/cipa_features.hpp>
#include <types/drug_block_input.hpp>
#include <types/mpi_profile.hpp>
#include <functions/inputoutput.hpp>

#include <cstring>
#include <cstdio>
#include <map>
#include <string>
#include <vector>

using std::multimap;
using std::string;
using std::vector;

int generate_report_drug(const Parameter *p_param, Cipa_Features *p_features, Drug_Block_Input &hill)
{
  
  char latex_filename[50];
  snprintf(latex_filename,sizeof(latex_filename),"report_drug.tex");
  FILE *fp_latex = fopen(latex_filename,"w");
  short font_size = 11;
  float image_scale = 0.8;
  const char *DOCUMENT_TITLE = "CardioSim Drug Simulation Result";
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
  fprintf(fp_latex,"Drug Name: %s\n\\\\Cmax: %s\n\\\\Ion Channel: 4 (INa, INaL, ICaL, IKr)\n", p_param->drug_name,p_param->concs);

  fprintf(fp_latex,"\\section*{Simulation Protocol}\n");
  fprintf(fp_latex,"Each drug was simulated by inducing %d beats at a %.0lf ms cycle.\n", p_param->pace_max, p_param->bcl);
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
  fprintf(fp_latex,"The biomarker for cardiotoxicity assessment proposed by the US FDA Leading Study Group is qNet, and Figure 3 shows the distribution of qNet for 1 compound (%s).\n", p_param->drug_name );
  fprintf(fp_latex,"\\\\The box-shaped range is the qNet data distribution interval between 5\\%% and 95\\%% confidence interval, the line interval is the standard deviation, and the dot distribution is the qNet data outside the standard deviation.\n");


  // get concentration values from input
  // and make the directories of them.
  std::vector<double> concs;
  char buffer[255];
  strncpy(buffer, p_param->concs, sizeof(buffer));
  char *token = strtok( buffer, "," );
  int idx = 0;
  while( token != NULL )
  { // begin data tokenizing
    concs.push_back(strtod(token, NULL));
        token = strtok(NULL, ",");
  } // end data tokenizing

  if(MPI_Profile::rank == 0) create_concs_directories(concs,p_param->drug_name);

  // generate time series for control
  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_Vm(mVolt).png}\n",
                    image_scale,p_param->drug_name, 0., p_param->drug_name, 0., 0 );
  fprintf(fp_latex,"\t\\caption{action potential result for control}\n");
  fprintf(fp_latex,"\t\\label{fig:image_vmplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, 0., 0);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_dVmdt(mVoltmsec).png}\n",
                    image_scale,p_param->drug_name, 0., p_param->drug_name, 0. , 0 );
  fprintf(fp_latex,"\t\\caption{rates of action potential result for control}\n");
  fprintf(fp_latex,"\t\\label{fig:image_dvmdtplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, 0., 0);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_Cai(x1000000)(nanoM).png}\n",
                    image_scale,p_param->drug_name, 0., p_param->drug_name, 0., 0 );
  fprintf(fp_latex,"\t\\caption{intracellular calcium concentration result for control}\n");
  fprintf(fp_latex,"\t\\label{fig:image_caiplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, 0., 0);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_ICaL(x1000)(nanoA).png}\n",
                    image_scale,p_param->drug_name, 0., p_param->drug_name, 0., 0 );
  fprintf(fp_latex,"\t\\caption{L-type calcium current result for control}\n");
  fprintf(fp_latex,"\t\\label{fig:image_icalplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, 0., 0);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_INaL(x1000)(nanoA).png}\n",
                    image_scale,p_param->drug_name, 0., p_param->drug_name, 0., 0 );
  fprintf(fp_latex,"\t\\caption{L-type sodium surrent result for control}\n");
  fprintf(fp_latex,"\t\\label{fig:image_inalplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, 0., 0);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_INa(x1000)(nanoA).png}\n",
                    image_scale,p_param->drug_name, 0., p_param->drug_name, 0.,0 );
  fprintf(fp_latex,"\t\\caption{sodium current result for control}");
  fprintf(fp_latex,"\t\\label{fig:image_inaplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, 0., 0);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_IKr(x1000)(nanoA).png}\n",
                    image_scale,p_param->drug_name, 0., p_param->drug_name, 0.,0 );
  fprintf(fp_latex,"\t\\caption{rapid potassium current result for control}\n");
  fprintf(fp_latex,"\t\\label{fig:image_icalplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, 0., 0);
  fprintf(fp_latex,"\\end{figure}\n");

  // attach time series result of the drug effect
  for( int sample_id = 0; sample_id < hill.size(); sample_id++ ){
    for( int jdx = 0; jdx < concs.size(); jdx++ ){
  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_Vm(mVolt).png}\n",
                    image_scale,p_param->drug_name, concs[jdx], p_param->drug_name, concs[jdx],sample_id );
  fprintf(fp_latex,"\t\\caption{action potential result for drug %s with concentration %.0lf of sample %d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\t\\label{fig:image_vmplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_dVmdt(mVoltmsec).png}\n",
                    image_scale,p_param->drug_name, concs[jdx], p_param->drug_name, concs[jdx],sample_id );
  fprintf(fp_latex,"\t\\caption{rates of action potential result for drug %s with concentration %.0lf of sample %d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\t\\label{fig:image_dvmdtplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_Cai(x1000000)(nanoM).png}\n",
                    image_scale,p_param->drug_name, concs[jdx], p_param->drug_name, concs[jdx],sample_id );
  fprintf(fp_latex,"\t\\caption{intracellular calcium concentration result for drug %s with concentration %.0lf of sample %d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\t\\label{fig:image_caiplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_ICaL(x1000)(nanoA).png}\n",
                    image_scale,p_param->drug_name, concs[jdx], p_param->drug_name, concs[jdx],sample_id );
  fprintf(fp_latex,"\t\\caption{L-type calcium current result for drug %s with concentration %.0lf of sample %d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\t\\label{fig:image_icalplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_INaL(x1000)(nanoA).png}\n",
                    image_scale,p_param->drug_name, concs[jdx], p_param->drug_name, concs[jdx],sample_id );
  fprintf(fp_latex,"\t\\caption{L-type sodium surrent result for drug %s with concentration %.0lf of sample %d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\t\\label{fig:image_inalplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_INa(x1000)(nanoA).png}\n",
                    image_scale,p_param->drug_name, concs[jdx], p_param->drug_name, concs[jdx],sample_id );
  fprintf(fp_latex,"\t\\caption{sodium current result for drug %s with concentration %.0lf of sample %d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\t\\label{fig:image_inaplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/%s_%.0lf_time_series_result_smp%d_plot_IKr(x1000)(nanoA).png}\n",
                    image_scale,p_param->drug_name, concs[jdx], p_param->drug_name, concs[jdx],sample_id );
  fprintf(fp_latex,"\t\\caption{rapid potassium current result for drug %s with concentration %.0lf of sample %d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\t\\label{fig:image_ikrplot_%s_%.0lf_%d}\n",
                    p_param->drug_name, concs[jdx], sample_id);
  fprintf(fp_latex,"\\end{figure}\n");

    } // end of concentration loop
  } // end of sample loop

  // generate figure of features distribution
  for(int jdx = 0; jdx < concs.size(); jdx++){
  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1lf]{./plots/%s/%.0lf/APD90.png}\n",
                    image_scale,p_param->drug_name, concs[jdx] );
  fprintf(fp_latex,"\t\\caption{APD90 distribution of %s with concentration %.0lf mMolar}\n",
                    p_param->drug_name, concs[jdx]);
  fprintf(fp_latex,"\t\\label{fig:image_apd90_%s_%.0lf}\n",
                    p_param->drug_name, 0.);
  fprintf(fp_latex,"\\end{figure}\n");

  fprintf(fp_latex,"\\begin{figure}[H]\n");
  fprintf(fp_latex,"\t\\centering\n");
  fprintf(fp_latex,"\t\\includegraphics[scale=%.1f]{./plots/%s/%.0lf/APD50.png}\n",
                    image_scale,p_param->drug_name, concs[jdx] );
  fprintf(fp_latex,"\t\\caption{APD50 distribution of %s with concentration %.0lf mMolar}\n",
                    p_param->drug_name, concs[jdx]);
  fprintf(fp_latex,"\t\\label{fig:image_apd90_%s_%.0lf}\n",
                    p_param->drug_name, 0.);
  fprintf(fp_latex,"\\end{figure}\n");
  }

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
