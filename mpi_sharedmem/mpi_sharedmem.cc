//Emanuel Gull 2016
#include<mpi.h>
#include<iostream>
#include<unistd.h>

int main(int argc, char **argv){

  MPI_Init(&argc, &argv);

  int global_rank, global_size;
  int host_rank, host_size;
  MPI_Comm_size(MPI_COMM_WORLD, &global_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);

  MPI_Comm hostcomm;
  MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &hostcomm);
  MPI_Comm_size(hostcomm, &host_size);
  MPI_Comm_rank(hostcomm, &host_rank);

  usleep(10000*global_rank);
  std::cout<<"global rank: "<<global_rank<<" global size: "<<global_size<<" host rank: "<<host_rank<<" host size: "<<host_size<<std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

  int memory_size=4096;
  MPI_Aint window_size=host_rank==0?memory_size:0;
  int disp_unit=sizeof(double);

  double *baseptr;

  MPI_Win win;
  MPI_Win_allocate_shared(window_size, disp_unit, MPI_INFO_NULL, hostcomm, &baseptr, &win);
  MPI_Win_shared_query(win, 0, &window_size, &disp_unit, &baseptr);

  MPI_Barrier(MPI_COMM_WORLD);
  if(host_rank==0) baseptr[0]=global_rank+42;
  MPI_Barrier(MPI_COMM_WORLD);

  usleep(10000*global_rank);
  std::cout<<"on global rank: "<<global_rank<<" local rank: "<<host_rank<<" read shared vector element: "<<baseptr[0]<<std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Win_free(&win);
  MPI_Finalize();
}
