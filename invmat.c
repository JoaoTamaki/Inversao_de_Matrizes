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

  printf("Sistema com a matriz original A e a matriz identidade B:\n"); // APAGAR
  prnSisLin(SL);          // APAGAR
  copia_matriz(SL->A, U, N);
  printf("U == A:\n");    // APAGAR
  prnMatriz(U, N);        // APAGAR

  // faz a fatoração LU e armazena as trocas de linha na LUT
  FatoracaoLU_PivoParcial(L, U, N, LUT, &tFatoracaoLU);
  printf("Tempo da fatoração LU: %lf\n\n", tFatoracaoLU);  // MUDAR PRA SAÍDA E APAGAR

  printf("U:\n");         // APAGAR
  prnMatriz(U, N);        // APAGAR
  printf("L:\n");         // APAGAR
  prnMatriz(L, N);        // APAGAR

  printf("LUT:\n");       // APAGAR
  prnVetorInt(LUT, N);    // APAGAR 

  // calcula o determinante para saber se a matriz é inversível (determinante != 0)
  real_t determinante = calculaDeterminante(U, N);
  printf("O determinante é %.15g.\n", determinante);        // APAGAR
  if ((determinante > -ERRO) && (determinante < ERRO)) {
    fprintf(stderr,"A matriz não é inversível.\n");
    exit(-3);
  }

  // cria copia
  SistLinear_t *SL_copia;
  SL_copia = alocaSisLin(SL->n, pontVet);
  copiaSisLin(SL, SL_copia);

  real_t **I;
  I = alocaMatriz(N);
  double tTotalY, tTotalX;
  int erro;
  erro = calculaInversa(L, U, I, LUT, N, &tTotalY, &tTotalX);

/*
  for (int i = 0; i < SL->n; i++){
    real_t *L_atual = (real_t*) malloc(N * N * sizeof(real_t)); 
    real_t *U1_atual = (real_t*) malloc(N * N * sizeof(real_t)); 

    real_t *Y_atual = (real_t*) malloc(N * N * sizeof(real_t)); 
    real_t *X_atual= (real_t*) malloc(N * N * sizeof(real_t)); 

    //Calcula Y:
    retrossubs(SL_copia, LUT, Y, i);
    CalculaYFROML(SL, LUT, Y);
    printf("Y:\n");
    prnMatriz(Y, SL_copia->n);

    //Calcula X:
    CalculaXFROMUY(X, U, Y, SL->n);
    printf("X:\n");
    prnMatriz(X, SL->n);

  }
    //libera tudo
    //free(X);
    //free(Y);
*/  

/*
  //Calcula L e U:
  fatoraLU(SL_copia, LUT, &tTotal);
  //Printa resultado: LU e Look up table
  printf("LU:\n");
  prnMatriz(SL_copia->A, SL->n);
  printf("LUT:\n");
  prnVetorInt (LUT, SL->n);

  //for (int iter = 0; iter < k; iter++){   
  if (k){
    //Alocação e teste
    real_t **Y, **X;
      
    Y = (real_t**) malloc(N * sizeof(real_t*));
    for(int i=0; i<N ; i++){
      Y[i] = (real_t*) malloc(N * sizeof(real_t));
    }
      
    X = (real_t**) malloc(N * sizeof(real_t*));
    for(int i=0; i<N ; i++){
      X[i] = (real_t*) malloc(N * sizeof(real_t));
    }
    //----------------------------------------INVERSAO DE MATRIZES----------------------------------------//
    for (int i = 0; i < SL->n; i++){

      //Calcula Y:
      //retrossubs(SL_copia, LUT, Y, i);
      //CalculaYFROML(SL, LUT, Y);
      //printf("Y:\n");
      //prnMatriz(Y, SL_copia->n);

      //Calcula X:
      //CalculaXFROMUY(X, U, Y, SL->n);
      //printf("X:\n");
      //prnMatriz(X, SL->n);

    }
    //libera tudo
    free(X);
    free(Y);
  }
  liberaSisLin(SL);
  liberaSisLin(SL_copia);
//n, t_egp, t_gs, it_gs, normaResiduo_gs, t_ref, it_ref, normaResiduo_ref
//SistLinear_t *SL, int *n, int tam_n, real_t *t_egp, real_t *t_gs, real_t *normaResiduo_gs, real_t *t_ref, int *it_ref, real_t *normaResiduo_ref
  if (!fp_in) fclose (fp_in);
  if (!fp_out) fclose (fp_out);
*/
  free(L);
  free(U);
  free(LUT);
}