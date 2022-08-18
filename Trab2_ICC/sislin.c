#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "sislin.h"

// alocaçao de matriz em memória 
SistLinear_t* alocaSisLin(unsigned int n) {
  SistLinear_t *SL = (SistLinear_t *) malloc(sizeof(SistLinear_t));
  
  if (SL) {    
    SL->n = n;
    SL->A = (real_t *) malloc(n * n * sizeof(real_t *));
    SL->b = (real_t *) malloc(n * n * sizeof(real_t *));

    if (!(SL->A) || !(SL->b)) {
      liberaSisLin(SL);
      return NULL;
    }
  }  
  return (SL);
}

// Liberacao de memória
void liberaSisLin(SistLinear_t *SL) {
  if (SL) {
    if (SL->A) {    
      free(SL->A);
    }
    // pequena alteração, pois b agora armazena a matriz identidade
    if (SL->b) {
      free(SL->b);
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
    
  if (tipo == hilbert) {
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j)  {
        SL->A[i*n+j] = 1.0 / (real_t)(i + j + 1);
      }
    }
  }
  else { // inicializa sistema normal e depois altera
    // inicializa a matriz A
    for (unsigned int i=0; i<n; ++i) {
      for (unsigned int j=0; j<n; ++j)  {
        SL->A[i*n+j] = (real_t)rand() * invRandMax;
      }
    }
    if (tipo == eqNula) {
      // sorteia eq a ser "nula"
      unsigned int nula = rand() % n;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[nula*n+j] = 0.0;
      }
    } 
    else if (tipo == eqProporcional) {
      // sorteia eq a ser "proporcional" e valor
      unsigned int propDst = rand() % n;
      unsigned int propSrc = (propDst + 1) % n;
      real_t mult = (real_t)rand() * invRandMax;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[propDst*n+j] = SL->A[propSrc*n+j] * mult;
      }
    } 
    else if (tipo == eqCombLinear) {
      // sorteia eq a ser "combLinear"
      unsigned int combDst = rand() % n;
      unsigned int combSrc1 = (combDst + 1) % n;
      unsigned int combSrc2 = (combDst + 2) % n;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[combDst*n+j] = SL->A[combSrc1*n+j] + SL->A[combSrc2*n+j];
      }
    }
    else if (tipo == diagDominante) {
      // aumenta o valor dos termos da diagonal principal
      for (unsigned int i=0; i<n; ++i) {
        real_t soma = 0.0;
        for (unsigned int j=0; j < i; ++j) soma += SL->A[i*n+j];
        for (unsigned int j=i+1; j < n; ++j) soma += SL->A[i*n+j];
        SL->A[i*n+i] += soma;
      }
    }
  }
}

SistLinear_t *lerSisLinArq(FILE *arqin) {
  unsigned int n;
  SistLinear_t *SL;
  fscanf(arqin, "%d",&n);

  SL = alocaSisLin(n);

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) 
      fscanf(arqin, "%lf", &SL->A[i*n+j]);
  
  return SL;
}

void prnSisLin(SistLinear_t *SL) {
  int n = SL->n;

  for(int i = 0; i < n; ++i) {
    printf("\n  ");
    for(int j = 0; j < n; ++j)
      printf("%.15g ", SL->A[i*n+j]);
  }
  printf("\n\n");
  for (int i = 0; i < n; ++i) {
    printf("\n  ");
    for (int j = 0; j < n; ++j)
      printf("%.15g ", SL->b[i*n+j]);
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

void prnMatriz(real_t *m, unsigned int n) {
  int i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; j++)
      printf ("%.15g \t", m[i*n+j]);
    printf ("\n");
  }
  printf ("\n\n");
}

void printaArquivoMatrizTransposta(FILE *fp_out, real_t *m, unsigned int n) {
  int i, j;

  for (i = 0; i < n; ++i) {
    for (j = 0; j < n; j++)
      fprintf (fp_out, "%.15g ", m[j*n+i]);
    fprintf (fp_out, "\n");
  }
  fprintf (fp_out, "\n\n");
}

void ordenaMatriz(real_t *m, real_t *mT, int *LUT, unsigned int n) {
  int i, j;
  for (i = 0; i < n; i ++){
    for (j = 0; j < n; j++)
      mT[LUT[i]*n+j] = m[i*n+j];
  }
}

real_t *alocaVetor(int n) {

  return (real_t *) malloc(n * sizeof(real_t));
}

int copia_matriz(real_t *x, real_t *y, int n) {

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) 
      y[i*n+j] = x[i*n+j];
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

      if (strcmp(argv[i],"-e") == 0){
        (*fp_in) = fopen(argv[i + 1], "r+");
      }else if (strcmp(argv[i],"-s") == 0){
        (*fp_out) = fopen(argv[i + 1], "w+");
      }else if (strcmp(argv[i],"-r") == 0){
        (*N) = atoi(argv[i + 1]);  
      }else if (strcmp(argv[i],"-i") == 0){
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

real_t* alocaMatriz(int N) {

  real_t *matriz;

  matriz = (real_t*) malloc(N * N * sizeof(real_t));
  if (!matriz) {
    fprintf(stderr,"Não foi possível alocar a matriz.\n");
    exit(-1);
  }
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