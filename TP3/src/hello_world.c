#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main(int argc, char** argv)
{
	char msg[20];
	int my_rank;
	int nbProc;
	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nbProc);

	if(my_rank %2){
		printf("Hello World ! My rank is %d. There are %d processes running. I'm odd, by the way.\n", my_rank, nbProc);
	} else {
		printf("Hi. I'm number %d. I'm even.\n", my_rank);
	}

	MPI_Finalize();
}

