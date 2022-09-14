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
  unsigned int flag_e, flag_s, flag_r, flag_i = 0;
  unsigned int N = 0;
  unsigned int k = 0;
  
  // verifica se o argumentos de entrada são válidos
  if (parseArguments(argc, argv, &fp_in, &fp_out, &N, &k, &flag_e, &flag_s, &flag_r, &flag_i) == -1)
    return 0;

  // aloca ponteiros e outras variáveis necessárias:
  SistLinear_t *SL;                         // SL, que vai armazenar a matriz original em A e a identidade em b
  real_t *L, *U;                            // matrizes L e U, que serão usadas na fatoração LU da matriz A
  int *LUT;                                 // look up table (LUT), que será usada para armazenar as trocas de linha no pivoteamento
  real_t tFatoracaoLU;                      // armazenará o tempo de execução da fatoração LU

  //inicializa o pad
  unsigned int pad;

  // faz a alocação do SL que vai armazenar a matriz conforme o tipo de entrada
  if (N > 0) {                              // caso a matriz seja gerada a partir de parâmetros
    SL = alocaSisLin(N, &pad);
    iniSisLin(SL, generico, COEF_MAX, pad);
  } else {                                  // caso a matriz seja dada via arquivo
    SL = lerSisLinArq(fp_in, &pad);
    N = SL->n;
  }
  criaMatrizIdentidade(SL->b, N, pad);      // inicializa os valores de b como uma matriz identidade
  L = alocaMatriz(N, pad);                  // aloca matriz L
  U = alocaMatriz(N, pad);                  // aloca matriz U
  LUT = alocaeInicilizaVetor(N, pad);       // aloca e inicializa LUT

  copia_matriz(SL->A, U, N, pad);

  // faz a fatoração LU e armazena as trocas de linha na LUT
  FatoracaoLU_PivoParcial(L, U, N, pad, LUT, &tFatoracaoLU);

  // calcula o determinante para saber se a matriz é inversível (determinante != 0)
  real_t determinante = calculaDeterminante(U, N, pad);
  if ((determinante > -ERRO) && (determinante < ERRO)) {
    fprintf(stderr,"A matriz não é inversível.\n");
    exit(-3);
  }

  //Inicializa variaveis para calculo da inversa e refinamento
  double op1, op2;                              //Tempo de cada iteração do calculo da iteração (Op1) e do resíduo (Op2)
  double tTempoResiduo = 0.0, tTempoIter = 0.0; //Tempo total do calculo da iteração (Op1) e do resíduo (Op2)
  real_t *I = alocaMatriz(N, pad);              //Matriz inversa conectada a LUT
  real_t *Inversa = alocaMatriz(N, pad);        //Matriz inversa ordenada para o refinamento
  real_t *x = alocaVetor(N, pad);               //Vetor X para o calculo do X
  real_t *y = alocaVetor(N, pad);               //Vetor Y para o calculo do Y
  real_t *norma = alocaVetor(k, pad);           //Vetor norma para armazenar as normas das iterações do calculo da inversa
  real_t *R = alocaMatriz(N, pad);              //Matriz resíduo
  real_t *W = alocaMatriz(N, pad);              //Matriz W para o refinamento
  real_t *vet = alocaVetor(N, pad);             //Vetor Auxiliar para o calculo de w

  //Calcula a primeira inversa que está linear em vez de colunar
  double tIter = timestamp();
  calculaInversa(L, U, I, x, y, LUT, N, pad);
  tIter = timestamp() - tIter;     
  //Ordena matriz inversa para ser usada sem LUT para o refinamento
  ordenaMatriz(I, Inversa, LUT, N, pad);

  //Realiza Refinamento k iterações
  for (int i = 0; i < k; i++){
    norma[i] = refinamento(SL, L, U, Inversa, R, W, vet, x, y, LUT, pad, &op1, &op2);
    tTempoIter += op1;
    tTempoResiduo += op2;
    fprintf(fp_out, "# iter %d: %.15g\n", i+1, norma[i]);
  }
  tTempoIter /= k;
  tTempoResiduo /= k;

  fprintf(fp_out, "# Tempo LU: %.15g\n", tFatoracaoLU);
  fprintf(fp_out, "# Tempo iter: %.15g\n", tTempoIter);
  fprintf(fp_out, "# Tempo residuo: %.15g\n", tTempoResiduo);
  fprintf(fp_out, "%d\n", N);
  printaArquivoMatrizTransposta(fp_out, Inversa, N, pad);

  //Libera ponteiros
  liberaSisLin(SL);
  free(L);
  free(U);
  free(LUT);
  free(I);
  free(x);
  free(y);
  free(Inversa);
  free(norma);
  free(R);
  free(W);
  free(vet);

  if (!fp_in) fclose (fp_in);
  if (!fp_out) fclose (fp_out);
}

//h35

// likwid-perfctr -C 3 -g L3 ./invmat -r 4 -i 10
// likwid-perfctr -C 3 -g L2CACHE ./invmat -r 4 -i 10
// likwid-perfctr -C 3 -g FLOPS_DP ./invmat -r 4 -i 10
// likwid-perfctr -C 3 -g FLOPS_AVX ./invmat -r 4 -i 10


/*! LISTA DE MELHORIAS
Duvidas:
  FEITO!

O que já era feito antes:
  -> Look up table
  -> Eliminar stride na multiplicação de matrizes 
  -> Questão de fazer matriz linear e colunar -> Fazer a tranposta apenas para a impressão para evitar Stride no acesso aos dados (multiplicação de matriz)

O que foi feito:
  -> Término do T1
  -> Matrizes alocadas contiguamente -> Melhorando o acesso a memória
  -> Economização de operações matemáticas
  -> Alligned_alloc -> Colunas não potências de 2 e linhas pode ser igual
  -> Padding -> Evitar cache trashing
  -> Restrict: Não aponta para o mesmo espaço de memória(vetor, matriz, ...)
  -> Inline: Tem que estar no mesmo codigo para expandir a função e eliminar o custo de empilhar parâmetros. 
  -> Instrumentar LIKWID nas duas versões
  -> Realizar loop unroll and jam

O que deve ser feito:
  -> Testar para onde usar restrict e inline (?)
  -> Script para testes

Observações:
  -> Por algum acaso no uroll and jam do calculo do x, fez ele variar um pouco (nas ultimas casas, ou seja, chega a ser irrelevante, mas não entendi porque ficou diferente)

Testes:
  -> N={32, 33, 64, 65, 128, 129, 256, 257, 512, 1000, 2000, 4000 6000 10000}
  -> -i 10
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
    ######L3 e L2CACHE
Apresentação:
  -> Realizar PDF com os gráficos e textos auxiliares
      O pacote deve ser arquivado e compactado com tar(1) e gzip(1) em um arquivo chamado login1-login2.tar.gz
      O arquivo contendo o relatório deve ser nomeado RELATORIO-login1-login2.pdf.
      Além disso, a conclusão deve apresentar os aspectos que sua equipe considerou mais relevantes/importantes no desenvolvimento do trabalho.
  -> No texto de entrega DEVE CONSTAR OBRIGATORIAMENTE o Nome e Números de Registro Acadêmico (RA) dos membros do grupo.
*/