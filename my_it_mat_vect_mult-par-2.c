/**
 *   \file my_it_mat_vect_mult.c
 *   \brief Multiplica iterativamente un matriz nxn 
 *          por un vector de n posiciones
 *
 *   \author Danny Múnera
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

const double MAXIMO = 25000;


/* función para generar <size> cantidad de datos aleatorios */
void gen_data(double * array, int size);
/* función para multiplicar iterativamente un matriz 
 * <m x n> por un vector de tam <n> */
void mat_vect_mult(double* A, double* x, double* y, int n, int it, int inicio, int final);
/* función para imprimir un vector llamado <name> de tamaño <m>*/
void print_vector(char* name, int rank, double*  z, int m);

int main()
{
  double* A = NULL;
  double* x = NULL;
  double* y = NULL;
  double* y_gather = NULL;
  int n=0, iters;
  long seed;
  int comm_sz;
  int rank;
  int inicio;
  int final;
  int segmentos;

  MPI_Init(NULL, NULL);

  // SPMD: rank , comm_sz
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("Número de procesos %d\n", comm_sz);

    if(rank==0){
        // Obtener las dimensiones
  printf("Ingrese la dimensión n:\n");
  scanf("%d", &n);
  printf("Ingrese el número de iteraciones:\n");
  scanf("%d", &iters);
  printf("Ingrese semilla para el generador de números aleatorios:\n");
  scanf("%ld", &seed);
  srand(seed);
  
}
  

  
  // Broadcast

  MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&iters, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&seed, 1, MPI_LONG, 0, MPI_COMM_WORLD);


  x = malloc(sizeof(double) * n);
  gen_data(x, n);
  A = malloc(sizeof(double) * n * n);
  gen_data(A, n*n);
  y = malloc(sizeof(double) * n);

    /*MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(A, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	  MPI_Barrier(MPI_COMM_WORLD);*/


    assert(n % comm_sz == 0);
    // local_n = n % comm_sz;
  // la matriz A tendrá una representación unidimensional
  

  

  //generar valores para las matrices
  
  
  inicio = rank*n;
  segmentos = (int)n/comm_sz;
  final=(rank+1)*n-1;
    
  printf("Proceso %d, inicio %d, final %d\n",rank, inicio, final);

  mat_vect_mult(A, x, y, n, iters, inicio, final);
  


  y_gather = malloc(sizeof(double) * n);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Allgather(&y, segmentos, MPI_DOUBLE, y_gather, n, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  print_vector("y", rank,y, n);

  free(A);
  free(x);
  free(y);
  


  MPI_Finalize();
  return 0;
}

void gen_data(double * array, int size){
  int i;
  for (i = 0; i < size; i++)
    array[i] = (double) rand() / (double) RAND_MAX;
}

void mat_vect_mult(double* A, double* x, double* y, int n, int it, int inicio, int fin){
  int h, i, j;
  for(h = 0; h < it; h++){  
    for(i = 0; i < n; i++){
      y[i] = 0.0;
      for(j = 0; j < n; j++){
	        y[i] += A[i*n+j] * x[j];
        }  
      }
    // x <= y
    for(i = 0; i < n; i++)
      x[i] = y[i];
  }
}

void print_vector(char* name,int rank, double*  z, int m) {
   int i;
   printf("\nProcess %d Vector %s\n",rank, name);
   for (i = 0; i < m; i++)
      printf("%f\n", z[i]);
   printf("\n");
}
