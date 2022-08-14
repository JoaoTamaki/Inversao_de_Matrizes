#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

int main (int argc, char** argv) {
  // inicializa gerador de números aleatórios
  srand(20221);

  // cria e inicializa variáveis de entrada do programa
  FILE *fp_out = stdout;
  FILE *fp_in = stdin;
  int flag_e, flag_s, flag_r, flag_i = 0;
  int N = 0;
  int k = 0;
  
  // verifica se o argumentos de entrada são válidos
  if (parseArguments(argc, argv, &fp_in, &fp_out, &N, &k, &flag_e, &flag_s, &flag_r, &flag_i) == -1)
    return 0;

  // aloca ponteiros e outras variáveis necessárias:
  SistLinear_t *SL;                     // SL, que vai armazenar a matriz original em A e a identidade em b
  real_t **L, **U;                      // matrizes L e U, que serão usadas na fatoração LU da matriz A
  int *LUT;                             // look up table (LUT), que será usada para armazenar as trocas de linha no pivoteamento
  real_t tFatoracaoLU;                  // armazenará o tempo de execução da fatoração LU

  // faz a alocação do SL que vai armazenar a matriz conforme o tipo de entrada
  if (N > 0) {                          // caso a matriz seja gerada a partir de parâmetros
    SL = alocaSisLin(N, pontVet);
    iniSisLin(SL, generico, COEF_MAX);
  } else {                              // caso a matriz seja dada via arquivo
    SL = lerSisLinArq(fp_in, generico);
    N = SL->n;
  }
  criaMatrizIdentidade(SL->b, N);       // inicializa os valores de b como uma matriz identidade
  L = alocaMatriz(N);                   // aloca matriz L
  U = alocaMatriz(N);                   // aloca matriz U
  LUT = alocaeInicilizaVetor(N);        // aloca e inicializa LUT

  copia_matriz(SL->A, U, N);

  // faz a fatoração LU e armazena as trocas de linha na LUT
  FatoracaoLU_PivoParcial(L, U, N, LUT, &tFatoracaoLU);

  // calcula o determinante para saber se a matriz é inversível (determinante != 0)
  real_t determinante = calculaDeterminante(U, N);
  if ((determinante > -ERRO) && (determinante < ERRO)) {
    fprintf(stderr,"A matriz não é inversível.\n");
    exit(-3);
  }

  printf("Matriz A:\n");
  printaMatriz(SL->A, N);

  real_t **I;
  I = alocaMatriz(N);
  double tTotalY, tTotalX;
  int erro;

  //Calcula a primeira inversa e arruma a ordem das linhas para ficar colunar para o calculo do resíduo.
  double tIter = timestamp();
  erro = calculaInversa(L, U, I, LUT, N, &tTotalY, &tTotalX);
  tIter += timestamp() - tIter;
  real_t **Inversa;
  Inversa = alocaMatriz(N);
  ordenaMatriz(I, Inversa, LUT, N);
  printf("Inversa\n");
  printaMatriz(Inversa, N);

  //Aloca vetor de norma
  real_t *norma;
  norma = alocaVetorZerado(norma, N);
  //Aloca matriz de resíduo
  real_t **R;                       //Matriz Resíduo R
  R = alocaMatriz(N);

  double tTempoResiduo = 0.0;
  for (int i = 0; i < k; i++){

    norma[i] = normaL2Residuo(SL, Inversa, R, N, &tIter);
    tTempoResiduo += tIter;
    fprintf(fp_out, "# iter %d: %.15g\n", i+1, norma[i]);
  }
  tTempoResiduo /= k;

  //prnVetor(norma, k);

  fprintf(fp_out, "# Tempo LU: %.15g\n", tFatoracaoLU);
  fprintf(fp_out, "# Tempo iter: %.15g\n", tIter);
  fprintf(fp_out, "Tempo residuo: %.15g\n", tTempoResiduo);
  fprintf(fp_out, "%d\n", N);
  prnMatriz(fp_out, I, N);

  liberaMatriz(I, N);
  liberaMatriz(Inversa, N);
  liberaMatriz(R, N);
  liberaMatriz(L, N);
  liberaMatriz(U, N);
  free(norma);
  free(LUT);
}