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

	MPI_Comm grid_comm;

	int dim_sizes[2];
	int wrap_around[2];
	int reorder = 1;
	dim_sizes[0] = 0;
	dim_sizes[1] = 0;
	wrap_around[0] = 1; 
	wrap_around[1] = 0;  
  
 	MPI_Dims_create(nbProc, 2, dim_sizes);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dim_sizes, wrap_around, reorder, &grid_comm);

	int coordinates[2]; 
	int my_grid_rank;
	MPI_Comm_rank ( grid_comm , &my_grid_rank ); 
	MPI_Cart_coords(grid_comm ,  my_grid_rank , 2, coordinates);

	printf("%d:[%d][%d]\n", my_rank, coordinates[0], coordinates[1]);

	MPI_Finalize();
	return EXIT_SUCCESS;
}

