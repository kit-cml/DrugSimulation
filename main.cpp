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
 
  // initialize MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &MPI_Profile::size );
  MPI_Comm_rank( MPI_COMM_WORLD, &MPI_Profile::rank );
  MPI_Get_processor_name( MPI_Profile::host_name, &MPI_Profile::host_name_len );

#if defined LICENSED
  int license_retval = 0;
  int check_license = 0;

  if (MPI_Profile::rank == 0) {
    check_license = validate_license();
    if (check_license < 0) {
      license_retval = -1;
    }
  }

  // Broadcast the license result to all processes
  MPI_Bcast(&license_retval, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // All processes check the license result
  if (license_retval < 0) {
    MPI_Finalize();
    return -1;
  }

  // Only proceed with barrier if license check passed
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  
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

  int err_code;
  double t_begin = MPI_Wtime();
  err_code = drug_simulation_bench(p_param,time_series,p_features,hill);
  if(err_code > 0) MPI_Abort(MPI_COMM_WORLD, err_code);
  double t_end = MPI_Wtime();

  // move to other post-processing app
  if(MPI_Profile::rank == 0) err_code = generate_report_drug(p_param);
  //MPI_Barrier(MPI_COMM_WORLD);

  if(err_code != 0) MPI_Abort(MPI_COMM_WORLD, err_code);
  mpi_printf(0,"Simulation finished at: %lf minutes.\n", (t_end-t_begin)*cml::math::SECONDS_TO_MINUTES);

  delete p_param;
  MPI_Finalize();

  return 0;
}
