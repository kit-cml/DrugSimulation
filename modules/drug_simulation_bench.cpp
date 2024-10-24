#include "drug_simulation_bench.hpp"
#include "insilico.hpp"
#include "postprocessing.hpp"

#include <types/cipa_features.hpp>
#include <types/drug_block_input.hpp>
#include <types/mpi_profile.hpp>
#include <functions/inputoutput.hpp>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iterator>
#include <map>
#include <mpi.h>
#include <string>
#include <vector>

using std::multimap;
using std::string;
using std::vector;

int drug_simulation_bench(const Parameter *p_param, multimap<double, string> &map_time_series, Cipa_Features &p_features, Drug_Block_Input &hill)
{
  char buffer[255];
  int error_code;

  // parsing ic50 and hill constants file into variables
  get_drug_data_from_file(p_param->hill_file, hill);
  error_code = check_drug_data_content(hill, p_param);
  if(error_code > 0) return 1;

  // parsing herg constants file into variables
  Drug_Block_Input herg;
  get_drug_data_from_file(p_param->herg_file, herg);
#ifdef ORD_DYN_2017
  error_code = check_drug_data_content(herg, p_param);
  if(error_code > 0) return 1;
  if(herg.size() != hill.size()) {
    mpi_printf(0, "Data size not equal. Hill: %d hERG %d\n", hill.size(), herg.size());
    return 1;
  }
#endif


  // get concentration values from input
  // and make the directories of them.
  std::vector<double> concs;
  concs.push_back(0.);
  strncpy(buffer, p_param->concs, sizeof(buffer));
  char *token = strtok( buffer, "," );
  while( token != NULL )
  { // begin data tokenizing
    concs.push_back(strtod(token, NULL));
    token = strtok(NULL, ",");
  } // end data tokenizing
  
  if(MPI_Profile::rank == 0) create_concs_directories(concs,p_param->drug_name);
  MPI_Barrier(MPI_COMM_WORLD);


  short sample_id = 0;
  int group_id = 0;
  bool is_ead;
  double inal_auc_control = 0.0;
  double ical_auc_control = 0.0;
  if( MPI_Profile::rank == 0 ){
    mpi_printf(0,"Running control in-silico simulation...\n");
    insilico(0., hill[sample_id], herg[sample_id], p_param, p_features, sample_id);
    mpi_printf(0,"Running control postprocessing simulation...\n");
    postprocessing(0., inal_auc_control, ical_auc_control, hill[sample_id], herg[sample_id], p_param, p_features, sample_id, group_id);
    inal_auc_control = p_features.inal_auc;
    ical_auc_control = p_features.ical_auc;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  printf("Before INaL at rank %d: %lf\n",MPI_Profile::rank, inal_auc_control);
  printf("Before ICaL at rank %d: %lf\n",MPI_Profile::rank, ical_auc_control);
  MPI_Bcast(&inal_auc_control, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&ical_auc_control, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  printf("After INaL at rank %d: %lf\n",MPI_Profile::rank, inal_auc_control);
  printf("After ICaL at rank %d: %lf\n",MPI_Profile::rank, ical_auc_control);

  for( sample_id = MPI_Profile::rank; sample_id < hill.size(); sample_id += MPI_Profile::size, group_id++ ){
    printf("Running sample# %d in Processor %d\n", sample_id, MPI_Profile::rank);
    for( short idx = 1; idx < concs.size(); idx++ ){
      insilico(concs[idx], hill[sample_id], herg[sample_id], p_param, p_features, sample_id);
      postprocessing(concs[idx], inal_auc_control, ical_auc_control, hill[sample_id], herg[sample_id], p_param, p_features, sample_id, group_id);
    }

  }
  MPI_Barrier(MPI_COMM_WORLD);

  return 0;
}
