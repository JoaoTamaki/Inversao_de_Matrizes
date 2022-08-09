#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

//----------------------------------------FUNCOES LAB1----------------------------------------// 

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
real_t normaL2Residuo(SistLinear_t *SL, real_t *I, int *LUT)
{
  real_t **R;
  R = alocaMatriz(SL->n);
  MultiplicaMQs(SL->A, I, R, LUT, SL->n);
  
  for (int i = 0; i < SL->n; i++){
    for (int j = 0; j < SL->n; j++){
      if (i == j){
        R[i][j] = 1.0 - R[i][j]; 
      }else{
        R[i][j] = 0.0 - R[i][j];
      }



    }
  }


  return sqrt(sum);
}

/*!
  \brief Método de Refinamento

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param erro menor erro aproximado para encerrar as iterações
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. Um nr positivo indica sucesso e o nr
          de iterações realizadas. Um nr. negativo indica um erro.
  */
/*int refinamento (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal, double *norma_ref)
{

  //*tTotal = timestamp();
  
  real_t *Ax = (real_t *) malloc(SL->n * sizeof(real_t));
  real_t *residuo = (real_t *) malloc(SL->n * sizeof(real_t));
  real_t *normaL2 = (real_t *) malloc(MAXIT * sizeof(real_t));
  real_t *save_b = (real_t *) malloc(SL->n * sizeof(real_t));

  copiaVetor(SL->b[0], save_b, SL->n);
  int numIter = 1;

  do{ //A*x=b -> A*A(t)*x=A(t)*b -> x=A(t)*b

    //calcula residuo (passo 2) para utilizar no passo 3
    for (int i = 0; i < SL->n; i++){
      for (int j = 0; j < SL->n; j++){
        Ax[i] = SL->A[i][j] * x[j];
      }
      residuo[i] = SL->b[i] - Ax[i];
    }

    copiaVetor(residuo, SL->b, SL->n);
    real_t *w = alocaVetorZerado(w, SL->n);
    real_t tIter;
    double norma_gs;
    int erro_w = eliminacaoGauss(SL, w, &tIter);
    copiaVetor(save_b, SL->b, SL->n);

    //calcula novo x (passo 4)
    for(int i = 0; i < SL->n; i++)
    {
      x[i] += w[i];
    }

    normaL2[numIter] = normaL2Residuo(SL, x);
    numIter++;
    free(w);
  }while(numIter < MAXIT && fabs(normaL2[numIter-1] - normaL2[numIter-2]) > ERRO);
  *norma_ref = normaL2[numIter-1];
  free(Ax);
  free(residuo);
  free(normaL2);

  //*tTotal = timestamp() - *tTotal;
  return (numIter-1);
}
*/

//----------------------------------------FUNÇÕES AUX----------------------------------------// 

/*!
  \brief Encontra o maior valor de uma coluna: pivô

  \param M ponteiro para a matriz
  \param pivNum índice da linha em análise que pode ser trocada pela linha pivô
  \param n tamanho da matriz

  \return índice da linha correspondente ao pivô (que tem o maior valor da coluna.
  */
int encontraMaxColunaPivo(double** M, int pivNum, int n) {

  int maxIndx = pivNum;
  double max = fabs(M[pivNum][pivNum]);

  for (int i = pivNum + 1; i < n; i++) {
    if (fabs(M[i][pivNum]) > max) {
      max = fabs(M[i][pivNum]);
      maxIndx = i;
    }
  }
  return maxIndx;
}

/*!
  \brief Inicializa valores de B como uma matriz identidade de mesma ordem que A

  \param M ponteiro para a matriz
  \param n tamanho da matriz

  \return B inicializado com a matriz identidade.
  */
void criaMatrizIdentidade(real_t **M, int n) {

  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if(i == j)
        M[i][j] = 1.0;
      else 
        M[i][j] = 0.0;
    }
  }

  return;
}

/*!
  \brief Esta função multiplica matrizes quadradas de mesma ordem. Ou seja, mR = mA x mB

  \param mA matriz da esquerda participante da multiplicação 
  \param mB matriz da direita participante da multiplicação
  \param mR matriz da resultante da multiplicação
  \param n  ordem das matrizes quadradas

  \return 

*/
void MultiplicaMQs(double** mA, double** mB, double** mR, int *LUT, int n) {

  for (int mri = 0; mri < n; mri ++) {
    for (int mrj = 0; mrj < n; mrj++) {
      mR[mri][mrj] = 0.0;
        for (int i = 0; i < n; i++)
          mR[mri][mrj] += mA[mri][i] * mB[LUT[i]][mrj];
    }
  }

  return;  
}

/*!
  \brief Calcula o determinante de uma matriz triangular superior

  \param M ponteiro para a matriz
  \param n tamanho da matriz

  \return o valor do determinante

*/
real_t calculaDeterminante(real_t **M, int n) {

  real_t determinante = 1.0;

  for (int i = 0;  i < n; ++i) 
    determinante = determinante * M[i][i];
  return determinante;
}

//----------------------------------------FUNÇÕES LU----------------------------------------// 
/*!
  \brief Fatoração LU. Decompõe a matriz U em L, U e P que guarda a ordem das trocas de linha.
  \Para isso, a função executa uma Eliminação de Gauss com pivoteamento parcial. 

  \param L matriz triangular inferior calculada com diagonal principal igual a 1
  \param U matriz triangular superior calculada
  \param P vetor que guarda as trocas de linhas durante o pivoteamento 
  \param n ordem das matrizes quadradas
  \param tTotal tempo total em milisegundos gastos pelo método

  \return código de erro. 0 em caso de sucesso.

*/

int FatoracaoLU_PivoParcial(double** L, double** U, int n, int* P, double *tTotal) {

  *tTotal = timestamp();

  // ESCALONAMENTO: PIVOTEAMENTO PARCIAL DA MATRIZ A COPIADA EM U, 
  // COM TRIANGULAÇÃO RESULTANDO EM MATRIZ SUPERIOR
  for (int piv = 0; piv < n; piv++) {               // para cada equação (linha) da matriz A 
    int newPiv = encontraMaxColunaPivo(U, piv, n);  // encontra o maior valor da coluna correspondente
    
    if (piv != newPiv) {                            // se o índice da linha pivô for diferente do índice da equação atual
      trocaLinhasMQ(U, n, piv, newPiv);             // troca linha de U atual pela linha pivô
      trocaLinhasMQ(L, n, piv, newPiv);             // troca linha de L atual pela linha pivô
      TrocaElementosVetor(P, piv, newPiv);          // armazena a troca na LUT
    }  
    L[piv][piv] = 1.0;                              // insere 1 na linha atual, na posição da diagonal principal de L
    double mPiv = U[piv][piv];

    for (int i = piv + 1; i < n; i++) {             // para cada coluna/ componente subsequente
      float mi = U[i][piv];
      U[i][piv] = 0.0;                              // vai zerando a parte superior da matriz U
      L[i][piv] = (mi/mPiv);                        // calcula o valor do multiplicador m
      if (isfinite(L[i][piv]) == 0) {
        fprintf(stderr, "O cálculo de um multiplicador de L na fatoração LU deu infinity ou NaN.\n");
        exit(-1);
      }
      for (int j = piv + 1; j < n; j++) {           // para o próximo elemento da linha da matriz U em diante
        U[i][j] -= L[i][piv] * U[piv][j];           // calcula o valor do novo coeficiente de U
        if (isfinite(U[i][j]) == 0) {
          fprintf(stderr, "O cálculo de um coeficiente de U na fatoração LU deu underflow ou overflow.\n");
          exit(-2);
        }
      }
    }
  }
    
  *tTotal = timestamp() - *tTotal;
  return 0;
}

void TrocaElementosVetor(int *vet, int i1, int i2){
  int aux = vet[i1];
  vet[i1] = vet[i2];
  vet[i2] = aux;
  return;
}

void trocaLinhasMQ(double **M, int n, int l1, int l2){
  double* aux;
  aux = M[l1];
  M[l1] = M[l2];
  M[l2] = aux;
  return;
}

/*!
  \brief Cácula Y da equação "L Y = I" como parte da fatoração LU.

  \param L matriz da fatoração LU. Matriz triangular inferior cuja diagonal principal consiste em somente 1's. 
  \param n tamanho n do sistema linear.
  \param LUT Look up table.
  \param n indice da Look up teble.
  \param t vetor y a ser calculado, que é usado para calcula Ux = y.

*/
void CalculaYFROML(real_t **L, int n, int *LUT, int k, real_t* y) {
  
  real_t *b = (real_t*)malloc(n * sizeof(real_t));
  for(int i = 0; i < n; i++) {
    b[i] = 0.0;
  }
  b[LUT[k]] = 1.0;
  
  for (int i = 0; i < n; i++) {
    y[i] = b[i];
    for (int j = 0; j < i; j++)
      y[i] -= L[i][j] * y[j];
    y[i] /= L[i][i];
  }

//  printf("L:\n");
//  prnMatriz(L, n);
  printf("y:\n");
  prnVetor(y, n);
//  printf("b:\n");
//  prnVetor(b, n);

  free(b);
  return;
}

/*!
  \brief Cácula X da equação "U X = Y" como parte da fatoração LU.

  \param U matriz da fatoração LU. Matriz triangular superior.  
  \param y vetor y previamente calculada na função CalculaYFROML.
  \param n  ordem das matrizes quadradas.
  \param x vetor x a ser calculado, que é uma linha da matriz identidade.
*/
void CalculaXFROMUY(real_t **U, real_t* y, int n, real_t *x){
  for (int i = n-1; i >= 0; --i) {
    x[i] = y[i];
    for (int j = i; j < n; ++j)
      x[i] -= U[i][j] * x[j];
    x[i] /= U[i][i];
  }
  return;
}

int calculaInversa(real_t **L, real_t **U, real_t **I, int *LUT, unsigned int n, double *tTotalY, double *tTotalX) {

  double tempo;
  int erro;

  real_t *y = (real_t*)malloc(n * sizeof(real_t));
  real_t *x = (real_t*)malloc(n * sizeof(real_t));
    
  for(int i = 0; i < n; i++) {

    tempo = timestamp();
    CalculaYFROML(L, n, LUT, i, y);
    *tTotalY += timestamp() - tempo;

    tempo = timestamp();
    CalculaXFROMUY(U, y, n, x);
    *tTotalX += timestamp() - tempo;

    for(int j = 0; j < n; j++) {
      I[i][j] = x[j];
    }
  }
  printf("I:\n");         // APAGAR
  prnMatriz(I, n);        // APAGAR

  free(y);
  free(x);

  *tTotalY /= n;
  *tTotalX /= n;
  return 0;
}