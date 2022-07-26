#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

int main (int argc, char** argv)
{
  // inicializa gerador de números aleatóreos
  srand(20221);

  //Inicializa as variaveis de entrada do programa
  FILE *fp_out = stdout;
  FILE *fp_in = stdin;
  int flag_e, flag_s, flag_r, flag_i = 0;
  int N = 0;
  int k = 0;
  
  if (parseArguments(argc, argv, &fp_in, &fp_out, &N, &k, &flag_e, &flag_s, &flag_r, &flag_i) == -1)
    return 0;

  //Aloca ponteiro para o sistema linear
  SistLinear_t *SL;

  //Tratamento da entrada
  if (N > 0){
    SL = alocaSisLin(N, pontVet);
    iniSisLin(SL, generico, COEF_MAX);
  }else{
    SL = lerSisLinArq (fp_in, generico);
  }

  //CRIA COPIA -> Acho que não precisa, pois temos a matriz L, U, Y e X da fatoração LU
  SistLinear_t *SL_copia;
  SL_copia = alocaSisLin(SL->n, pontVet);
  copiaSisLin(SL, SL_copia);

  //INICIALIZA B
  MatrizIndentidade(SL->b, SL->n);
  
  printf("A:\n");
  prnMatriz(SL->A, SL->n);

  printf("B:\n");
  prnMatriz(SL->b, SL->n);

  //Create Look up table
  int *LUT = (int*) malloc (SL->n * sizeof(int));
  for (int i = 0; i < SL->n; i++){
    LUT[i] = i;
  }

  real_t tTotal; //Inicializa o tempo para cada iteração -> Não sei aonde começa o calculo
  //Calcula L e U:
  fatoraLU(SL_copia, LUT, &tTotal);
  //Printa resultado: LU e Look up table
  printf("LU:\n");
  prnMatriz(SL_copia->A, SL->n);
  printf("LUT:\n");
  prnVetorInt (LUT, SL->n);


  for (int iter = 0; iter < k; iter++){   
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
}


