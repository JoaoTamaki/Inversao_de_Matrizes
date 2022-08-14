#ifndef __SISLIN_H__
#define __SISLIN_H__

#define COEF_MAX 32.0 // Valor máximo usado para gerar valores aleatórios de
		      // coeficientes nos sistemas lineares.

// Tipo de alocação para matrizes
typedef enum {
  pontPont=0, // Matriz como vetor de N ponteiros para vetores de tamanho N
  pontVet     // Matriz como vetor de N ponteiros para um único vetor de tamanho N*N
} tipoAloc_t;

// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  real_t **A; // coeficientes
  real_t **b; // termos independentes
  unsigned int n; // tamanho do SL
  tipoAloc_t tipoAloc_A; // tipo de alocação usada na matriz de coeficientes
} SistLinear_t;

// Tipos de matrizes de coeficientes usados pela função 'inicializaSistLinear()'
typedef enum {
    generico = 0,
    hilbert,
    diagDominante,
    eqNula,
    eqProporcional,
    eqCombLinear
} tipoSistLinear_t;


// Alocaçao e desalocação de matrizes
SistLinear_t* alocaSisLin (unsigned int n, tipoAloc_t tipo);
real_t *alocaVetorZerado(real_t *x, int n);
void liberaSisLin (SistLinear_t *SL);
void liberaMatriz(real_t **m, unsigned int n);
void iniSisLin (SistLinear_t *SL, tipoSistLinear_t tipo, real_t coef_max);

// Leitura e impressão de sistemas lineares
SistLinear_t *lerSisLinArq (FILE *arqin, tipoAloc_t tipo);
void prnSisLin (SistLinear_t *SL);
void prnVetorInt (int *v, unsigned int n);
void prnVetor (real_t *v, unsigned int n);
void prnMatriz(FILE *fp_out, real_t **m, unsigned int n);
void printaMatriz(real_t **m, unsigned int n);
void ordenaMatriz(real_t **m, real_t **mT, int *LUT, unsigned int n);

// Outras funções úteis
int copia_matriz(real_t **x, real_t **y, int n);
int copia_vetor(real_t *x, real_t *y, int n);
int copiaSisLin(SistLinear_t *SL, SistLinear_t *SL_copia);
int parseArguments(int argc, char** argv, FILE** fp_in, FILE** fp_out, int *N, int *k, int *flag_e, int *flag_s, int *flag_r, int *flag_i);
real_t** alocaMatriz(int N);
int* alocaeInicilizaVetor(int N);

#endif // __SISLIN_H__