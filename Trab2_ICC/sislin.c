#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "sislin.h"

// Alocaçao e desalocação de estruturas
SistLinear_t* alocaSisLin(unsigned int n, unsigned int *pad) {
  //Definição do PAD
  if (n % 2 == 0)
    *pad = 1;
  else
    *pad = 2;

  printf("Pad: %d\n", *pad);

  SistLinear_t *SL = aligned_alloc(16, (2*n*(*pad))*sizeof(real_t) + sizeof(SistLinear_t));    //Ta certo?
  if (SL) {    
    SL->n = n;
    SL->A = alocaMatriz(n, *pad);
    SL->b = alocaMatriz(n, *pad);
    if (!(SL->A) || !(SL->b)) {
      liberaSisLin(SL);
      return NULL;
    }
  }  
  return (SL);
}

real_t* alocaVetor(unsigned int n, unsigned int pad) {

  real_t *vetor;
  vetor = aligned_alloc(16, (n+pad)*sizeof(real_t));
  if (!vetor) {
    fprintf(stderr,"Não foi possível alocar o vetor.\n");
    exit(-1);
  }
  return (vetor);
}

int* alocaVetorInt(unsigned int n, unsigned int pad) {

  int *vetor;
  unsigned int np = n + pad;
  vetor = aligned_alloc(16, (np)*sizeof(int));
  if (!vetor) {
    fprintf(stderr,"Não foi possível alocar o vetor.\n");
    exit(-1);
  }
  return (vetor);
}

real_t* alocaMatriz(unsigned int n, unsigned int pad) {

  real_t *matriz;
  unsigned int np = n + pad;
  matriz = aligned_alloc(16, n*np*sizeof(real_t));
  if (!matriz) {
    fprintf(stderr,"Não foi possível alocar a matriz.\n");
    exit(-1);
  }
  return (matriz);
}

int* alocaeInicilizaVetor(unsigned int n, unsigned int pad) {

  int *vetor;
  vetor = alocaVetorInt(n, pad);
  for (int i = 0; i < n; i++)
    vetor[i] = i;

  return (vetor);
}

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
void iniSisLin(SistLinear_t *SL, tipoSistLinear_t tipo, real_t coef_max, unsigned int pad) {
  unsigned int n = SL->n;
  unsigned int np = n+pad;
  // para gerar valores no intervalo [0,coef_max]
  real_t invRandMax = ((real_t)coef_max / (real_t)RAND_MAX);
    
  if (tipo == hilbert) {
    for (unsigned int i = 0; i < n; ++i) {
      for (unsigned int j = 0; j < n; ++j)  {
        SL->A[i*np+j] = 1.0 / (real_t)(i + j + 1);
      }
    }
  }
  else { // inicializa sistema normal e depois altera
    // inicializa a matriz A
    for (unsigned int i=0; i<n; ++i) {
      for (unsigned int j=0; j<n; ++j)  {
        SL->A[i*np+j] = (real_t)rand() * invRandMax;
      }
    }
    if (tipo == eqNula) {
      // sorteia eq a ser "nula"
      unsigned int nula = rand() % n;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[nula*np+j] = 0.0;
      }
    } 
    else if (tipo == eqProporcional) {
      // sorteia eq a ser "proporcional" e valor
      unsigned int propDst = rand() % n;
      unsigned int propSrc = (propDst + 1) % n;
      real_t mult = (real_t)rand() * invRandMax;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[propDst*np+j] = SL->A[propSrc*np+j] * mult;
      }
    } 
    else if (tipo == eqCombLinear) {
      // sorteia eq a ser "combLinear"
      unsigned int combDst = rand() % n;
      unsigned int combSrc1 = (combDst + 1) % n;
      unsigned int combSrc2 = (combDst + 2) % n;
      for (unsigned int j=0; j<n; ++j) {
        SL->A[combDst*np+j] = SL->A[combSrc1*np+j] + SL->A[combSrc2*np+j];
      }
    }
    else if (tipo == diagDominante) {
      // aumenta o valor dos termos da diagonal principal
      for (unsigned int i=0; i<n; ++i) {
        real_t soma = 0.0;
        for (unsigned int j=0; j < i; ++j) soma += SL->A[i*np+j];
        for (unsigned int j=i+1; j < n; ++j) soma += SL->A[i*np+j];
        SL->A[i*np+i] += soma;
      }
    }
  }
}

//Leitura e escrita de estruturas

SistLinear_t *lerSisLinArq(FILE *arqin, unsigned int *pad) {
  unsigned int n;
  SistLinear_t *SL;
  fscanf(arqin, "%d",&n);
  unsigned int np = n+(*pad);
  SL = alocaSisLin(n, pad);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) 
      fscanf(arqin, "%lf", &SL->A[i*np+j]);
  
  return SL;
}

void prnSisLin(SistLinear_t *SL, unsigned int pad) {
  unsigned int n = SL->n;
  unsigned int np = n + pad;

  prnMatriz(SL->A, n, pad);
  prnMatriz(SL->b, n, pad);
}

void prnVetorInt(int *v, unsigned int n) {

  printf ("\n");
  for (int i = 0; i < n; ++i)
    printf ("%d ", v[i]);
  printf ("\n\n");
}

void prnVetor(real_t *v, unsigned int n) {

  printf ("\n");
  for(int i = 0; i < n; ++i)
      printf ("%.15g ", v[i]);
  printf ("\n\n");
}

void prnMatriz(real_t *m, unsigned int n, unsigned int pad) {

  unsigned int np = n + pad;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; j++)
      printf ("%.15g \t", m[i*np+j]);
    printf ("\n");
  }
  printf ("\n\n");
}

void printaArquivoMatrizTransposta(FILE *fp_out, real_t *m, unsigned int n, unsigned int pad) {

  unsigned int np = n + pad;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; j++)
      fprintf (fp_out, "%.15g ", m[j*np+i]);
    fprintf (fp_out, "\n");
  }
  fprintf (fp_out, "\n\n");
}

//Outras funções uteis

void ordenaMatriz(real_t *m, real_t *mT, int *LUT, unsigned int n, unsigned int pad) {

  unsigned int np = n + pad;
  for (int i = 0; i < n; i ++){
    for (int j = 0; j < n; j++)
      mT[LUT[i]*np+j] = m[i*np+j];
  }
}

int copia_matriz(real_t *x, real_t *y, unsigned int n, unsigned int pad) {
  unsigned int np = n + pad;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) 
      y[i*np+j] = x[i*np+j];
  }
  return 0;
}

int copia_vetor(real_t *x, real_t *y, unsigned int n, unsigned int pad) {

  for(int i = 0; i < n; i++)
    y[i] = x[i];

  return 0;
}

int copiaSisLin(SistLinear_t *SL, SistLinear_t *SL_copia, unsigned int pad) {

  SL_copia->n = SL->n;
  copia_matriz(SL->A, SL_copia->A, SL->n, pad);
  copia_matriz(SL->b, SL_copia->b, SL->n, pad);
  return 0;
}

int parseArguments(int argc, char** argv, FILE** fp_in, FILE** fp_out, unsigned int *N, unsigned int *k, unsigned int *flag_e, unsigned int *flag_s, unsigned int *flag_r, unsigned int *flag_i) {

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