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

	uint32_t tabSend[N];
	if(my_rank == 0){
		for(int i = 0; i < N; i++){
			tabSend[i] = ((((int)tabSend / N) * (i+1)) % (N*10))/8;
		}
	} else {
		for(int i = 0; i < N; i++){
			tabSend[i] = (((int)tabSend / N) * (i+1)) %4; // initializes to 0
 		}
	}	
    MPI_Bcast(tabSend, N, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

	for(int i = 0; i < N; i++){
        printf("process %d : tab[%d] = %u\n", my_rank, i, tabSend[i]);
    }
    printf("\n");

	MPI_Finalize();
	return EXIT_SUCCESS;
}

