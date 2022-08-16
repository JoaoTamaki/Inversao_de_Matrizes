#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

/*! LISTA DE MELHORIAS
Duvidas:
  -> O que deve ser feito de mudança para cada iteração do resíduo?
  -> Como usar alligned_alloc
  -> Questão de fazer matriz linear e colunar -> Stride
  -> Revisar restrict e inline
  -> Qual a necessidade de fazer padding?
  -> No texto de entrega DEVE CONSTAR OBRIGATORIAMENTE o Nome e Números de Registro Acadêmico (RA) dos membros do grupo. (???)
  
O que já era feito antes:
  -> Look up table

Gerais:
  -> Matrizes alocadas contiguamente -> Melhorando o acesso a memória
  -> Economização de operações matemáticas

Possíveis melhorias:
  -> void *aligned_alloc(alignment, size)
  -> Verificar stride
  -> Verificar restrict
  -> Verificar inline
  -> Padding
  -> Realizar loop unroll and jam

O que deve ser feito:
  -> Tirar dúvidas
  -> Instrumentar LIKWID nas duas versões
  -> Realizar loop unroll and jam
  -> N={32, 33, 64, 65, 128, 129, 256, 257, 512, 1000, 2000, 4000 6000 10000}
  -> -i 10
  -> Matrizes sempre aleatórias
  -> Ver como realizar os gráficos
      Cada teste deve apresentar em linhas distintas os valores para o cálculo de cada operações (op1 e op2). Assim, os gráficos terão sempre 4 linhas, duas para a v1 e duas para a v2;
      Cada gráfico deve ser explicado e você deve demonstrar que consegue entender o que está reportado nele;
      Os gráficos devem ser apresentados com o eixo das ordenadas (eixo y) em escala logarítmica.

      Os seguintes testes devem ser executados (um gráfico para cada teste):
        Teste de tempo: mostra o tempo médio do cálculo da op1 e o tempo médio do cálculo da op2 (utilize a função "timestamp()" para medir o tempo);
        Banda de Memória: utilizar o grupo MEM do LIKWID, e apresentar o resultado de "Memory bandwidth [MBytes/s]"; Caso não tenha o grupo MEM, utilize o grupo L3.
        Cache miss L1: utilizar o grupo CACHE ou L1CACHE do LIKWID, e apresentar o resultado de "data cache miss ratio".
        Caso não tenha o cache miss da L1, utilize o cache miss da L2 (grupo L2CACHE)
        Operações aritméticas: utilizar o grupo FLOPS_DP ou FLOPS_AVX do LIKWID, e apresentar o resultado de "MFLOP/s"
  -> Realizar PDF com os gráficos e textos auxiliares
      O pacote deve ser arquivado e compactado com tar(1) e gzip(1) em um arquivo chamado login1-login2.tar.gz
      O arquivo contendo o relatório deve ser nomeado RELATORIO-login1-login2.pdf.
      Além disso, a conclusão deve apresentar os aspectos que sua equipe considerou mais relevantes/importantes no desenvolvimento do trabalho.

*/



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
  real_t *L, *U;                        // matrizes L e U, que serão usadas na fatoração LU da matriz A
  int *LUT;                             // look up table (LUT), que será usada para armazenar as trocas de linha no pivoteamento
  real_t tFatoracaoLU;                  // armazenará o tempo de execução da fatoração LU

  // faz a alocação do SL que vai armazenar a matriz conforme o tipo de entrada
  if (N > 0) {                          // caso a matriz seja gerada a partir de parâmetros
    SL = alocaSisLin(N);
    iniSisLin(SL, generico, COEF_MAX);
  } else {                              // caso a matriz seja dada via arquivo
    SL = lerSisLinArq(fp_in);
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
  prnMatriz(SL->A, N);

  real_t *I;
  I = alocaMatriz(N);
  double tTotalY, tTotalX;
  int erro;

  //Calcula a primeira inversa e arruma a ordem das linhas para ficar colunar para o calculo do resíduo.
  double tIter = timestamp();
  erro = calculaInversa(L, U, I, LUT, N, &tTotalY, &tTotalX);           //OP 2
  tIter += timestamp() - tIter;
  real_t *Inversa;
  Inversa = alocaMatriz(N);
  ordenaMatriz(I, Inversa, LUT, N);
  printf("Inversa\n");
  prnMatriz(Inversa, N);

  //Aloca vetor de norma
  real_t *norma;
  norma = alocaVetorZerado(norma, N);
  //Aloca matriz de resíduo
  real_t *R;                       //Matriz Resíduo R
  R = alocaMatriz(N);

  double tTempoResiduo = 0.0;
  for (int i = 0; i < k; i++){

    norma[i] = normaL2Residuo(SL, Inversa, R, N, &tIter);               //OP 2
    tTempoResiduo += tIter;
    fprintf(fp_out, "# iter %d: %.15g\n", i+1, norma[i]);
  }
  tTempoResiduo /= k;

  //prnVetor(norma, k);

  fprintf(fp_out, "# Tempo LU: %.15g\n", tFatoracaoLU);
  fprintf(fp_out, "# Tempo iter: %.15g\n", tIter);
  fprintf(fp_out, "Tempo residuo: %.15g\n", tTempoResiduo);
  fprintf(fp_out, "%d\n", N);
  printaArquivoMatrizTransposta(fp_out, Inversa, N);

  free(I);
  free(Inversa);
  free(R);
  free(L);
  free(U);
  free(SL);
  free(norma);
  free(LUT);
}