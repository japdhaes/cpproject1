#include "mainapplication.h"
#include <mpi.h>

int main(int argc, char *argv[])
{
//    argc=1;
//    argv=new char*[0];
    int numprocs, myrank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MainApplication m(myrank, numprocs);
    m.runApplication();

    MPI_Finalize();
    return 0;
}

