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

/* función para generar <size> cantidad de datos aleatorios */
void gen_data(double * array, int size);
/* función para multiplicar iterativamente un matriz 
 * <m x n> por un vector de tam <n> */
void mat_vect_mult(double* A, double* x, double* y, int n, int it);
/* función para imprimir un vector llamado <name> de tamaño <m>*/
void print_vector(char* name, double*  y, int m);

int main()
{
  double* A = NULL;
  double* x = NULL;
  double* y = NULL;
  int n, iters;
  long seed;
  int comm_sz;
  int rank;

    double local_n;
  MPI_Init(NULL, NULL);

  // SPMD: rank , comm_sz
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  if(rank == 0){
    // Obtener las dimensiones
  printf("Ingrese la dimensión n:\n");
  scanf("%d", &n);
  printf("Ingrese el número de iteraciones:\n");
  scanf("%d", &iters);
  printf("Ingrese semilla para el generador de números aleatorios:\n");
  scanf("%ld", &seed);
  srand(seed);
  x = malloc(sizeof(double) * n);
  gen_data(x, n);
  }
    assert(n % comm_sz == 0);
    local_n = n % comm_sz;
    MPI_Bcast(&x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);


  // la matriz A tendrá una representación unidimensional
  A = malloc(sizeof(double) * local_n * n);
  y = malloc(sizeof(double) * local_n);

  //generar valores para las matrices
  gen_data(A, local_n*n);
  
  mat_vect_mult(A, x, y, local_n, iters);

  print_vector("y", y, local_n);
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

void mat_vect_mult(double* A, double* x, double* y, int n, int it){
  int h, i, j;
  for(h = 0; h < it; h++){
    for(i = 0; i < n; i++){
      y[i] = 0.0;
      for(j = 0; j < n; j++)
	y[i] += A[i*n+j] * x[j];
    }
    // x <= y
    for(i = 0; i < n; i++)
      x[i] = y[i];
  }
}

void print_vector(char* name, double*  y, int m) {
   int i;
   printf("\nVector %s\n", name);
   for (i = 0; i < m; i++)
      printf("%f ", y[i]);
   printf("\n");
}
