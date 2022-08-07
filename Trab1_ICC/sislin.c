#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "sislin.h"

// alocaçao de matriz em memória 
SistLinear_t* alocaSisLin(unsigned int n, tipoAloc_t tipo) {
  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
  
  if (SL) {    
    SL->n = n;
    SL->tipoAloc_A = tipo;
    SL->A = (real_t **) malloc(n * sizeof(real_t *));

    // alteração na alocação de B para armazenar a matriz identidade
    SL->b = (real_t **) malloc(n * sizeof(real_t *));
    for (int i = 0; i < n; i++)
      SL->b[i] = (real_t *) malloc(n * sizeof(real_t));

    if (!(SL->A) || !(SL->b)) {
      liberaSisLin(SL);
      return NULL;
    }

    // Matriz como vetor de N ponteiros para um único vetor com N*N elementos
    if (tipo == pontVet) {
      SL->A[0] = (real_t *) malloc(n * n * sizeof(real_t));
      if (!(SL->A[0])) {
        liberaSisLin(SL);
        return NULL;
      }	
      for (int i = 1; i < n; ++i)
        SL->A[i] = SL->A[i-1]+n;
    }
    else if (tipo == pontPont) {  // Matriz  como  vetor de  N  ponteiros
				                          // para N vetores de N elementos cada
      for (int i = 0; i < n; ++i)
        SL->A[i] = (real_t *) malloc(n * sizeof(real_t));
    }
  }
  return (SL);
}

// Liberacao de memória
void liberaSisLin(SistLinear_t *SL) {
  if (SL) {
    if (SL->A) {
      if (SL->tipoAloc_A == pontVet) {
        if (SL->A[0]) free (SL->A[0]);
      }
      else if (SL->tipoAloc_A == pontPont) {
        for (int i = 0; i < SL->n; ++i) free (SL->A[i]);
      }      
      free(SL->A);
    }
    
    // pequena alteração, pois b agora armazena a matriz identidade
    if (SL->b[0]) {
      free(SL->b[0]);
      for (int i = 0; i < SL->n; ++i) free (SL->A[i]);
    }
    free(SL);
  }
}

/*!
  \brief Cria coeficientes e termos independentes do SL. FOI COMENTADO O B, JÁ QUE NO TRABALHO SERÁ A MATRIZ IDENTIDADE
  *
  \param SL Ponteiro para o sistema linear
  \param tipo Tipo de sistema linear a ser criado. Pode ser: 
     comSolucao, eqNula, eqProporcional, eqCombLinear, hilbert 
  \param coef_max Maior valor para coeficientes e termos independentes
*/
void iniSisLin(SistLinear_t *SL, tipoSistLinear_t tipo, real_t coef_max) {
  unsigned int n = SL->n;
  // para gerar valores no intervalo [0,coef_max]
  real_t invRandMax = ((real_t)coef_max / (real_t)RAND_MAX);

  // inicializa vetor b
  //for (unsigned int i = 0; i < n; ++i) {
    //SL->b[i] = (real_t)rand() * invRandMax;
  //}
    
  if (tipo == hilbert) {
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j)  {
        SL->A[i][j] = 1.0 / (real_t)(i + j + 1);
      }
    }
  }
  else { // inicializa sistema normal e depois altera
    // inicializa a matriz A
    for (unsigned int i=0; i<n; ++i) {
      for (unsigned int j=0; j<n; ++j)  {
        SL->A[i][j] = (real_t)rand() * invRandMax;
      }
    }
    if (tipo == eqNula) {
      // sorteia eq a ser "nula"
      unsigned int nula = rand() % n;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[nula][j] = 0.0;
      }
    } 
    else if (tipo == eqProporcional) {
      // sorteia eq a ser "proporcional" e valor
      unsigned int propDst = rand() % n;
      unsigned int propSrc = (propDst + 1) % n;
      real_t mult = (real_t)rand() * invRandMax;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[propDst][j] = SL->A[propSrc][j] * mult;
      }
    } 
    else if (tipo == eqCombLinear) {
      // sorteia eq a ser "combLinear"
      unsigned int combDst = rand() % n;
      unsigned int combSrc1 = (combDst + 1) % n;
      unsigned int combSrc2 = (combDst + 2) % n;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[combDst][j] = SL->A[combSrc1][j] + SL->A[combSrc2][j];
      }
    }
    else if (tipo == diagDominante) {
      // aumenta o valor dos termos da diagonal principal
      for (unsigned int i=0; i<n; ++i) {
        real_t soma = 0.0;
        for (unsigned int j=0; j < i; ++j) soma += SL->A[i][j];
        for (unsigned int j=i+1; j < n; ++j) soma += SL->A[i][j];
        SL->A[i][i] += soma;
      }
    }
  }
}

SistLinear_t *lerSisLinArq(FILE *arqin, tipoAloc_t tipo) {
  unsigned int n;
  SistLinear_t *SL;
  fscanf(arqin, "%d",&n);

  SL = alocaSisLin(n, tipo);

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) 
      fscanf(arqin, "%lf", &SL->A[i][j]);
  
  return SL;
}

void prnSisLin(SistLinear_t *SL) {
  int n=SL->n;

  for(int i = 0; i < n; ++i) {
    printf("\n  ");
    for(int j = 0; j < n; ++j)
      printf("%.15g ", SL->A[i][j]);
  }
  printf("\n\n");
  for (int i = 0; i < n; ++i) {
    printf("\n  ");
    for (int j = 0; j < n; ++j)
      printf("%.15g ", SL->b[i][j]);
  }
  printf("\n\n");
}

void prnVetorInt(int *v, unsigned int n) {
  int i;

  printf ("\n");
  for (i = 0; i < n; ++i)
    printf ("%d ", v[i]);
  printf ("\n\n");
}

void prnVetor(real_t *v, unsigned int n) {
  int i;

  printf ("\n");
  for(i = 0; i < n; ++i)
      printf ("%.15g ", v[i]);
  printf ("\n\n");
}

void prnMatriz(real_t **m, unsigned int n) {
  int i, j;

  printf ("\n");
  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; j++)
      printf ("%.15g ", m[i][j]);
    printf ("\n");
  }
  printf ("\n\n");
}

real_t *alocaVetorZerado(real_t *x, int n) {

  x = (real_t *) malloc(n * sizeof(real_t));
  for (int i = 0; i < n; i++)
    x[i] = 0.0;
  return x;
}

int copia_matriz(real_t **x, real_t **y, int n) {

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) 
      y[i][j] = x[i][j];
  }
  return 0;
}

int copia_vetor(real_t *x, real_t *y, int n) {

  for(int i = 0; i < n; i++)
    y[i] = x[i];

  return 0;
}

int copiaSisLin(SistLinear_t *SL, SistLinear_t *SL_copia) {

  SL_copia->n = SL->n;
  copia_matriz(SL->A, SL_copia->A, SL->n);
  copia_matriz(SL->b, SL_copia->b, SL->n);
  return 0;
}

int parseArguments(int argc, char** argv, FILE** fp_in, FILE** fp_out, int *N, int *k, int *flag_e, int *flag_s, int *flag_r, int *flag_i) {

  if (argc >= 3 || argc <= 9) {
    for (int i = 1; i < argc; i = i + 2) {

      if (strcmp(argv[i],"-e") == 0)
        (*fp_in) = fopen(argv[i + 1], "r+");
      else if (strcmp(argv[i],"-s") == 0)
        (*fp_out) = fopen(argv[i + 1], "w+");
      else if (strcmp(argv[i],"-r") == 0)
        (*N) = atoi(argv[i + 1]);  
      else if (strcmp(argv[i],"-i") == 0) {
        (*k) = atoi(argv[i + 1]);
        //validação da condição: k > 0
        if ((*k) <= 0){
          fprintf(stderr,"Entrada inválida: %s\n", argv[i]);
          return -1;
        }      
      } else {
        fprintf(stderr,"Argumento inválido %s\n", argv[i]);
        return -1;
      }
    }
    return 0;
  }
  fprintf(stderr,"Número inválido de argumentos. Confira as entradas possíveis no README.\n");
  return -1;
}

real_t** alocaMatriz(int N) {

  real_t **matriz;

  matriz = (real_t**) malloc(N * sizeof(real_t*));
  if (!matriz) {
    fprintf(stderr,"Não foi possível alocar a matriz.\n");
    exit(-1);
  }
  for (int i = 0; i < N; i++)
    matriz[i] = (real_t*) malloc(N * sizeof(real_t));

  return (matriz);
}


int* alocaeInicilizaVetor(int N) {

  int *vetor;

  vetor = (int*) malloc (N * sizeof(int));
  if (!vetor) {
    fprintf(stderr,"Não foi possível alocar o vetor.\n");
    exit(-1);
  }
  for (int i = 0; i < N; i++)
    vetor[i] = i;

  return (vetor);

}