#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

//----------------------------------------FUNÇÕES AUX----------------------------------------//
void calculaNovoI(real_t *I, real_t *W, unsigned int n, unsigned int pad){
  int i, j;
  unsigned int mult, np = n + pad;
  for (i = 0; i < n; i++){
    mult = i*np;
    for (j = 0; j < n; j++){
      I[mult+j] += W[mult+j]; 
    }
  }
}

real_t normaL2Residuo(real_t *R, unsigned int n, unsigned int pad){
  unsigned int i, j;
  unsigned int mult, np = n + pad;
  real_t sum;
  for (int i = 0; i < n; i++){
    mult = i*np;               //em vez de executar i*j vezes, faz apenas i vezes
    for (int j = 0; j < n; j++){
      sum += R[mult+j]*R[mult+j];
    }
  }
  return sqrt(sum);
}

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 
  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
real_t refinamento(SistLinear_t *SL, real_t *L, real_t *U, real_t *I, real_t *R, real_t *W, real_t *vet, real_t *x, real_t *y, int *LUT, unsigned int pad, double *tTempoIter, double *tTempoResiduo)
{
  //Inicializa variáveis
  unsigned int n = SL->n;
  unsigned int np = n + pad;
  real_t normaResiduo;
  
  *tTempoIter = timestamp();  
  //Realiza Resíduo = B - A * I -> B(Identidade), A(Matriz de entrada), I(Matriz Inversa)
  calculaMatrizResiduo(SL->A, SL->b, I, R, n, pad);
  //Aw = R
  calculaW(L, U, R, W, vet, x, y, LUT, n, pad);
  calculaNovoI(I, W, n, pad);
  *tTempoIter = timestamp() - *tTempoIter;

  //Calcula norma do resíduo
  *tTempoResiduo = timestamp();
  normaResiduo = normaL2Residuo(R, n, pad);
  *tTempoResiduo = timestamp() - *tTempoResiduo;
  return normaResiduo;
}

/*!
  \brief Esta função multiplica matrizes quadradas de mesma ordem. Ou seja, mR = mA x mB
  \param mA   matriz A participante da multiplicação 
  \param mB   matriz B participante da subtração
  \param mI   matriz Inversa participante da multiplicação
  \param mR   matriz Resíduo resultante
  \param n    ordem das matrizes quadradas
  \param pad  tamanho do pad
  \return 
*/
void calculaMatrizResiduo(real_t* mA, real_t* mB, real_t* mI, real_t* mR, unsigned int n, unsigned int pad) {

  unsigned int mult, np = n + pad;
  for (int i = 0; i < n; i ++) {
    mult = i*np;               //em vez de executar i*j vezes, faz apenas i vezes
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < n; k++)
        mR[mult+j] += mA[mult+k] * mI[j*np+k];
    }
    mR[mult+i] = mB[mult+i] - mR[mult+i]; 
  }
  return;  
}

/*!
  \brief Encontra o maior valor de uma coluna: pivô
  \param M      ponteiro para a matriz
  \param pivNum índice da linha em análise que pode ser trocada pela linha pivô
  \param n      tamanho da matriz
  \param pad    tamanho do pad
  \return       índice da linha correspondente ao pivô (que tem o maior valor da coluna)
  */
unsigned int encontraMaxColunaPivo(double* M, unsigned int pivNum, unsigned int n, unsigned int pad) {

  int i;
  unsigned int np = n + pad;
  unsigned int maxIndx = pivNum;
  double max = fabs(M[pivNum*np+pivNum]);
  unsigned int mult;                               //em vez de executar 2*i vezes, faz apenas i vezes caso entre no (IF)

  for (i = pivNum + 1; i < n; i++) {
    mult = i*np+pivNum;
    if (fabs(M[mult]) > max) {
      max = fabs(M[mult]);
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
void criaMatrizIdentidade(real_t *M, unsigned int n, unsigned int pad) {

  int i, j;
  unsigned int np = n + pad;

  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if(i == j)
        M[i*np+j] = 1.0;
      else 
        M[i*np+j] = 0.0;
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
real_t calculaDeterminante(real_t *M, unsigned int n, unsigned int pad) {

  real_t determinante = 1.0;
  int i;
  unsigned int np = n + pad;
  for (i = 0;  i < n; ++i) 
    determinante = determinante * M[i*np+i];
  return determinante;
}

//----------------------------------------FUNÇÕES LU----------------------------------------// 
/*!
  \brief Fatoração LU. Decompõe a matriz U em L, U e P que guarda a ordem das trocas de linha. Para isso, a função executa uma Eliminação de Gauss com pivoteamento parcial. 
  \param L matriz triangular inferior calculada com diagonal principal igual a 1
  \param U matriz triangular superior calculada
  \param LUT vetor que guarda as trocas de linhas durante o pivoteamento 
  \param n ordem das matrizes quadradas
  \param tTotal tempo total em milisegundos gastos pelo método
  \return código de erro. 0 em caso de sucesso.
*/
int FatoracaoLU_PivoParcial(double* L, double* U, unsigned int n, unsigned int pad, int* LUT, double *tTotal) {

  *tTotal = timestamp();
  int i, j;
  unsigned int piv, newPiv, np = n + pad;
  real_t mPiv, mi;
  // ESCALONAMENTO: PIVOTEAMENTO PARCIAL DA MATRIZ A COPIADA EM U, 
  // COM TRIANGULAÇÃO RESULTANDO EM MATRIZ SUPERIOR
  for (piv = 0; piv < n; piv++) {                       // para cada equação (linha) da matriz A 
    newPiv = encontraMaxColunaPivo(U, piv, n, pad);     // encontra o maior valor da coluna correspondente
    if (piv != newPiv) {                                // se o índice da linha pivô for diferente do índice da equação atual
      trocaLinhasMQ(U, n, pad, piv, newPiv);            // troca linha de U atual pela linha pivô
      trocaLinhasMQ(L, n, pad, piv, newPiv);            // troca linha de L atual pela linha pivô
      TrocaElementosVetor(LUT, piv, newPiv);            // armazena a troca na LUT
    } 
    L[piv*np+piv] = 1.0;                                 // insere 1 na linha atual, na posição da diagonal principal de L
    mPiv = U[piv*np+piv];

    for (i = piv + 1; i < n; i++) {                     // para cada coluna/ componente subsequente
      mi = U[i*np+piv];
      U[i*np+piv] = 0.0;                                // vai zerando a parte superior da matriz U
      L[i*np+piv] = (mi/mPiv);                          // calcula o valor do multiplicador m
      if (isfinite(L[i*np+piv]) == 0) {
        fprintf(stderr, "O cálculo de um multiplicador de L na fatoração LU deu infinity ou NaN.\n");
        exit(-1);
      }
      for (j = piv + 1; j < n; j++) {                   // para o próximo elemento da linha da matriz U em diante
        U[i*np+j] -= L[i*np+piv] * U[piv*np+j];         // calcula o valor do novo coeficiente de U
        if (isfinite(U[i*np+j]) == 0) {
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

void trocaLinhasMQ(double *M, unsigned int n, unsigned int pad, int l1, int l2){
  double aux;
  int j;
  unsigned int mult1, mult2, np = n + pad; 
  for (j = 0; j < n; j++){
    //em vez de 2*n operações, realiza n vezes
    mult1 = l1*np+j;
    mult2 = l2*np+j;

    aux = M[mult1];
    M[mult1] = M[mult2];
    M[mult2] = aux;
  }
  return;
}

/*!
  \brief Calcula Y da equação "L Y = I" como parte da fatoração LU.
  \param L matriz da fatoração LU. Matriz triangular inferior cuja diagonal principal consiste em somente 1's. 
  \param n tamanho n do sistema linear.
  \param LUT Look up table.
  \param k índice da Look up table.
  \param y vetor y a ser calculado, que é usado para calcula Ux = y.
*/
void CalculaYFROML(real_t *L, unsigned int n, unsigned int pad, int *LUT, int k, real_t* y, real_t *b) {
  int i, j;
  unsigned int mult, np = n + pad;
  for (i = 0; i < n; i++) {               //percorre diagonal principal de cima pra baixo
    y[i] = b[i];                          //inicia o y com o b da linha
    mult = i*np;                          //em vez de executar i*j vezes, faz apenas i vezes
    for (j = 0; j < i; j++)               //j inicia no inicio e percorre antes da diagonal principal
      y[i] -= L[mult+j] * y[j];           //realiza subtrações de Lx de cada posição após a diagonal principal y = b[i] - L[i][i+1]*y[i] - ... - L[i][n]*y[i]
  }
  return;
}

/*!
  \brief Calcula X da equação "U X = Y" como parte da fatoração LU.
  \param U matriz da fatoração LU. Matriz triangular superior.  
  \param y vetor y previamente calculada na função CalculaYFROML.
  \param n ordem das matrizes quadradas.
  \param x vetor x a ser calculado, que é uma linha da matriz identidade.
*/
void CalculaXFROMUY(real_t *U, real_t* y, unsigned int n, unsigned int pad, real_t *x){
  int i, j;
  unsigned int mult, np = n + pad;
  for (i = n-1; i >= 0; --i) {            //percorre diagonal principal de baixo pra cima
    x[i] = y[i];                          //inicia o x com o y da linha
    mult = i*np;                          //em vez de executar i*(np-j+1)+1 vezes, faz apenas i+1 vezes
    for (j = i+1; j < n; ++j){            //j inicia após a diagonal principal 
      x[i] -= U[mult+j] * x[j];           //realiza subtrações de Uy de cada posição após a diagonal principal x = y[i] - U[i][i+1]*x[i] - ... - U[i][n]*x[i]
    }
    x[i] /= U[mult+i];                    //divide pelo elemento da diagonal principal de U[i][i]
  }
  return;
}

/*!
  \brief Calcula matriz inversa usando fatoração LU e eliminação de Gauss 
  \param L matriz da fatoração LU. Matriz triangular inferior cuja diagonal principal consiste em somente 1's.
  \param U matriz da fatoração LU. Matriz triangular superior.  
  \param I matriz que armazenará a inversa de A.
  \param LUT Look up table que armazena as trocas de linha em b.
  \param n tamanho n do sistema linear.
  \param t vetor y a ser calculado, que é usado para calcula Ux = y.
*/
void calculaInversa(real_t *L, real_t *U, real_t *I, real_t *x, real_t *y, int *LUT, unsigned int n, unsigned int pad) {
  int i, j;
  unsigned int np = n + pad;
  real_t *b = aligned_alloc(16, (np)*sizeof(real_t)); 
  for(i = 0; i < n; i++) {
    //inicializa b
    for(j = 0; j < n; j++) {    
      b[j] = 0.0;
    }
    b[LUT[i]] = 1.0;
    CalculaYFROML(L, n, pad, LUT, i, y, b);
    CalculaXFROMUY(U, y, n, pad, x);
    for(j = 0; j < n; j++) {
      I[LUT[i]*np+j] = x[j];
    }
  }
  free(b);
}

void calculaW(real_t *L, real_t *U, real_t *R, real_t *W, real_t *vet, real_t *x, real_t *y, int *LUT, unsigned int n, unsigned int pad) {
  int i, j;
  unsigned int mult, np = n + pad;
  for(i = 0; i < n; i++) {
    mult = i*np;                          //em vez de executar i*j vezes, faz apenas i vezes
    //Invete R para usar no LUw=R
    for (j = 0; j < n; j++){
      vet[LUT[j]] = R[mult+j];
    }
    CalculaYFROML(L, n, pad, LUT, i, y, vet);
    CalculaXFROMUY(U, y, n, pad, x);
    mult = LUT[i]*np;
    for(j = 0; j < n; j++) {
      W[mult+j] = x[j];
    }
  }
}