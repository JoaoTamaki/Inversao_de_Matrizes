#ifndef __SISLIN_H__
#define __SISLIN_H__

#define COEF_MAX 32.0 // Valor máximo usado para gerar valores aleatórios de
		                  // coeficientes nos sistemas lineares.

// Tipo de alocação para matrizes
typedef enum {
  pontPont=0,     // Matriz como vetor de N ponteiros para vetores de tamanho N
  pontVet         // Matriz como vetor de N ponteiros para um único vetor de tamanho N*N
} tipoAloc_t;

// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  real_t *A;      // coeficientes
  real_t *b;      // termos independentes
  unsigned int n; // tamanho do SL
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


// Alocaçao e desalocação de estruturas
SistLinear_t* alocaSisLin(unsigned int n, unsigned int *pad);
real_t *alocaVetor(unsigned int n, unsigned int pad);
int* alocaVetorInt(unsigned int n, unsigned int pad);
real_t* alocaMatriz(unsigned int n, unsigned int pad);
int* alocaeInicilizaVetor(unsigned int n, unsigned int pad);
void liberaSisLin(SistLinear_t *SL);
void iniSisLin(SistLinear_t *SL, tipoSistLinear_t tipo, real_t coef_max, unsigned int pad);

// Leitura e impressão de sistemas lineares
SistLinear_t *lerSisLinArq(FILE *arqin, unsigned int *pad);
void prnSisLin(SistLinear_t *SL, unsigned int pad);
void prnVetorInt (int *v, unsigned int n);
void prnVetor (real_t *v, unsigned int n);
void prnMatriz(real_t *m, unsigned int n, unsigned int pad);
void printaArquivoMatrizTransposta(FILE *fp_out, real_t *m, unsigned int n, unsigned int pad);

// Outras funções úteis
void ordenaMatriz(real_t *m, real_t *mT, int *LUT, unsigned int n, unsigned int pad);
int copia_matriz(real_t *x, real_t *y, unsigned int n, unsigned int pad);
int copia_vetor(real_t *x, real_t *y, unsigned int n, unsigned int pad);
int copiaSisLin(SistLinear_t *SL, SistLinear_t *SL_copia, unsigned int pad);
int parseArguments(int argc, char** argv, FILE** fp_in, FILE** fp_out, unsigned int *N, unsigned int *k, unsigned int *flag_e, unsigned int *flag_s, unsigned int *flag_r, unsigned int *flag_i);

#endif // __SISLIN_H__