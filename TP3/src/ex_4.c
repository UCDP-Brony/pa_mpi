#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <stdint.h>
#include <stdlib.h>

int main(int argc, char** argv)
{
	int my_rank;
	int nbProc;
	MPI_Status status;

	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nbProc);

	if(nbProc != 2){
		printf("attention, il faut utiliser deux processus. Pas plus, pas moins.\n");
		return EXIT_FAILURE;
	}

	if(my_rank == 0){
		uint64_t tabSend[10];
		for(int i = 0; i < 10; i++){
			tabSend[i] = 0xffffffffffffffff / (10 - i);
		}
		MPI_Send(tabSend,10,MPI_UNSIGNED_LONG, 1, 0xbeef, MPI_COMM_WORLD);
		MPI_Recv(tabSend,10,MPI_UNSIGNED_LONG, 1, 0xbeef, MPI_COMM_WORLD, &status);
		
	} else if(my_rank == 1){
		uint64_t tabRecept[10];
		MPI_Recv(tabRecept,10,MPI_UNSIGNED_LONG, 0, 0xbeef, MPI_COMM_WORLD, &status);
		for(int i = 0; i < 10; i++){
			tabRecept[i]-=1;
		}
		MPI_Send(tabRecept,10,MPI_UNSIGNED_LONG, 0, 0xbeef, MPI_COMM_WORLD);	
	}
	MPI_Finalize();
	return EXIT_SUCCESS;
}

