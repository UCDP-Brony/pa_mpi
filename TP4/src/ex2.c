#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdint.h>
#include <stdlib.h>

#define N 10

int main(int argc, char** argv)
{
	int my_rank;
	int nbProc;
	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nbProc);

    uint32_t sum = 0;
    MPI_Reduce(&my_rank, &sum, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&sum, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	printf("process %d : our sum is %d\n", my_rank, sum);

	MPI_Finalize();
	return EXIT_SUCCESS;
}

