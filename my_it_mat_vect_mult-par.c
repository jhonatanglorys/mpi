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
void print_vector(char* name, int rank, double*  z, int m);

int main()
{
  double* A = NULL;
  double* x = NULL;
  double* y = NULL;
  double* fila=NULL;
  double* subtotal=NULL;
  int datos;
  int n, iters;
  long seed;
  int comm_sz;
  int rank;



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

    // la matriz A tendrá una representación unidimensional
    A = malloc(sizeof(double) * n * n);
    x = malloc(sizeof(double) * n);
    y = malloc(sizeof(double) * n);

    //generar valores para las matrices
    gen_data(A, n*n);
    gen_data(x, n);
  }

    
    assert(n % comm_sz == 0);

    fila = malloc(sizeof(double) * n);
    subtotal = malloc(sizeof(double) * n);
    
    datos = (int)n/comm_sz;

    MPI_Scatter(A, datos, MPI_DOUBLE, fila, datos, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(x,n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(y,n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD);
    mat_vect_mult(fila, x, subtotal, n, iters);

    MPI_Gather(subtotal, 1, MPI_DOUBLE, y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   
    print_vector("fila", rank,fila, datos);
   
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

void print_vector(char* name,int rank, double*  z, int m) {
   int i;
   printf("\nProcess %d Vector %s\n",rank, name);
   for (i = 0; i < m; i++)
      printf("%f\n", z[i]);
   printf("\n");
}
