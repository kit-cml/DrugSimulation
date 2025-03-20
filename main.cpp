#include <functions/inputoutput.hpp>
#include <functions/softdongle.hpp>
#include <types/cipa_features.hpp>
#include <types/cml_consts.hpp>
#include <types/mpi_profile.hpp>
#include <types/parameter.hpp>
#include "modules/drug_simulation_bench.hpp"
#include "modules/report_drug.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mpi.h>

using std::multimap;

int main(int argc, char **argv)
{
  // enable real-time output in stdout
  setvbuf( stdout, NULL, _IONBF, 0 );
  // Set the locale to one that uses a comma as the decimal separator
  //setlocale(LC_NUMERIC, "de_DE"); // German locale (or any other locale that uses a comma)
  multimap <double, string> time_series;
 
#if defined LICENSED
  int check_license = validate_license();
  if(check_license < 0){
    return -1;
  }
#endif
  // initialize MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &MPI_Profile::size );
  MPI_Comm_rank( MPI_COMM_WORLD, &MPI_Profile::rank );
  MPI_Get_processor_name( MPI_Profile::host_name, &MPI_Profile::host_name_len );

  // input parameter object
  Parameter *p_param;
  p_param = new Parameter();
  p_param->init();
  assign_params(&argc,argv,p_param);
  p_param->show_val();

  // cipa feature object to store
  // the 14 features and time series 
  // to the report
  Cipa_Features p_features;

  // drug input data
  Drug_Block_Input hill;

  double t_begin = MPI_Wtime();
  drug_simulation_bench(p_param,time_series,p_features,hill);
  double t_end = MPI_Wtime();

  // move to other post-processing app
  if(MPI_Profile::rank == 0) generate_report_drug(p_param);
  // MPI_Barrier(MPI_COMM_WORLD);

  mpi_printf(0,"Simulation finished at: %lf minutes.\n", (t_end-t_begin)*cml::math::SECONDS_TO_MINUTES);

  delete p_param;
  MPI_Finalize();

  return 0;
}
