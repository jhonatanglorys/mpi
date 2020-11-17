#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char** argv){
	// MPI
	int PE, nPE;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nPE);
	MPI_Comm_rank(MPI_COMM_WORLD, &PE);
	if(PE==0){
		MPI_Status status;
		int node;
		for(int i=1; i<nPE; i++){
			MPI_Recv(&node, 1, MPI_INT, i, 10, MPI_COMM_WORLD, &status);
		}
	}else{
		MPI_Send(&PE, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	/*
	 * Only one process, i.e. process 0, will execute this command.
	 */
	if(PE == 0){
		printf("Using %d processes.\n", nPE);
	}

	/*
	 * All processes will exectute this command, although each have a different PE
	 * number.
	 */
	printf("Process %d reporting.\n", PE);

	// Parameters
	int numElements = 10;
	int *elements;
	elements = (int*) calloc(numElements, sizeof(double));

	// Calculation
    	for(int i=0+PE; i<numElements; i=i+nPE){
		elements[i] = i;
	}

	// Output
	MPI_Barrier(MPI_COMM_WORLD);
    	for(int i=0; i<nPE; i++){
		if(PE == i){
			printf("Process %d copy of elements pre-broadcast\n", PE);
    		for(int j=0; j<numElements; j++){
			printf("%d ", elements[j]);
		}
		printf("\n");
		}
	}

	// Broadcast
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i=0; i<numElements; i++){
		MPI_Bcast(&elements[i], 1, MPI_INT, i%nPE, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

    	for(int i=0; i<nPE; i++){
		if(PE == i){
			printf("Process %d copy of elements post-broadcast\n", PE);
    		for(int j=0; j<numElements; j++){
			printf("%d ", elements[j]);
		}
		printf("\n");
		}
	}

	// Closing
	free(elements);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
    return 0;
}