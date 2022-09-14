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

  //Inicializa variaveis para calculo da inversa e refinamento
  double op1, op2;                              //Tempo no calculo de x e y
  double tTempoResiduo = 0.0, tTempoIter = 0.0; //Tempo total do calculo da iteração (Op1) e do resíduo (Op2)
  real_t **I = alocaMatriz(N);                   //Matriz inversa conectada a LUT
  real_t **Inversa = alocaMatriz(N);             //Matriz inversa ordenada para o refinamento
  real_t *x = alocaVetor(N);                    //Vetor X para o calculo do X
  real_t *y = alocaVetor(N);                    //Vetor Y para o calculo do Y
  real_t *norma = alocaVetor(k);                //Vetor norma para armazenar as normas das iterações do calculo da inversa
  real_t **R = alocaMatriz(N);                   //Matriz resíduo
  real_t **W = alocaMatriz(N);                   //Matriz W para o refinamento
  real_t *vet = alocaVetor(N);                  //Vetor Auxiliar para o calculo de w

  //Calcula a primeira inversa e arruma a ordem das linhas para ficar colunar para o calculo do resíduo.
  double tIter = timestamp();
  calculaInversa(L, U, I, x, y, LUT, N);
  tIter = timestamp() - tIter;
  //Ordena matriz inversa para ser usada sem LUT para o refinamento
  ordenaMatriz(I, Inversa, LUT, N);

  //Realiza Refinamento k iterações
  for (int i = 0; i < k; i++){
    norma[i] = refinamento(SL, L, U, Inversa, R, W, vet, x, y, LUT, &op1, &op2);
    tTempoIter += op1;
    tTempoResiduo += op2;
    fprintf(fp_out, "# iter %d: %.15g\n", i+1, norma[i]);
  }
  tTempoIter /= k;
  tTempoResiduo /= k;

  fprintf(fp_out, "# Tempo LU: %.15g\n", tFatoracaoLU);
  fprintf(fp_out, "# Tempo iter: %.15g\n", tTempoIter);
  fprintf(fp_out, "Tempo residuo: %.15g\n", tTempoResiduo);
  fprintf(fp_out, "%d\n", N);
  printaArquivoMatrizTransposta(fp_out, Inversa, N);

  //Libera ponteiros
  liberaSisLin(SL);
  liberaMatriz(L, N);
  liberaMatriz(U, N);
  free(LUT);
  liberaMatriz(I, N);
  free(x);
  free(y);
  liberaMatriz(Inversa, N);
  free(norma);
  liberaMatriz(R, N);
  liberaMatriz(W, N);
  free(vet);

  if (!fp_in) fclose (fp_in);
  if (!fp_out) fclose (fp_out);
}