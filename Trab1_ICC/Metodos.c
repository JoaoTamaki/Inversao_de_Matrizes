#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>

#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

//----------------------------------------FUNÇÕES AUX----------------------------------------// 
/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 
  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
void calculaNovoI(real_t **I, real_t **W, unsigned int n){
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      I[i][j] += W[i][j]; 
    }
  }
}

real_t normaL2Residuo(real_t **R, unsigned int n){
  
  int sum;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      sum += R[i][j]*R[i][j];
    }
  }
  return sqrt(sum);
}

real_t refinamento(SistLinear_t *SL, real_t **L, real_t **U, real_t **I, real_t **R, real_t **W, real_t *vet, real_t *x, real_t *y, int *LUT, double *tTempoIter, double *tTempoResiduo){
  
  //Inicializa variáveis
  unsigned int n = SL->n;
  real_t normaResiduo;

  *tTempoIter = timestamp();  
  //Realiza Resíduo = B - A * I -> B(Identidade), A(Matriz de entrada), I(Matriz Inversa)
  calculaMatrizResiduo(SL->A, SL->b, I, R, n);
  calculaW(L, U, R, W, vet, x, y, LUT, n);
  calculaNovoI(I, W, n);
  *tTempoIter += timestamp() - *tTempoIter;

  //Calcula norma do resíduo
  *tTempoResiduo = timestamp();
  normaResiduo = normaL2Residuo(R, n);
  *tTempoResiduo += timestamp() - *tTempoResiduo;

  return normaResiduo;
}

/*!
  \brief Esta função multiplica matrizes quadradas de mesma ordem. Ou seja, mR = mA x mB
  \param mA matriz da esquerda participante da multiplicação 
  \param mB matriz da direita participante da multiplicação
  \param mR matriz da resultante da multiplicação
  \param n  ordem das matrizes quadradas
  \return 
*/
void calculaMatrizResiduo(real_t** mA, real_t** mB, real_t** mI, real_t** mR, int n) {

  for (int i = 0; i < n; i ++) {
    for (int j = 0; j < n; j++) {
      mR[i][j] = 0.0;
      for (int k = 0; k < n; k++)
        mR[i][j] += mA[i][k] * mI[j][k];
    }
    mR[i][i] = mB[i][i] - mR[i][i]; 
  }
  return;  
}

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
  \brief Calcula Y da equação "L Y = I" como parte da fatoração LU.
  \param L matriz da fatoração LU. Matriz triangular inferior cuja diagonal principal consiste em somente 1's. 
  \param n tamanho n do sistema linear.
  \param LUT Look up table.
  \param k índice da Look up table.
  \param y vetor y a ser calculado, que é usado para calcula Ux = y.
*/
void CalculaYFROML(real_t **L, int n, int *LUT, int k, real_t* y, real_t *b) {
  for (int i = 0; i < n; i++) {     //percorre diagonal principal de cima pra baixo
    y[i] = b[i];                    //inicia o y com o b da linha
    for (int j = 0; j < i; j++)     //j inicia no inicio e percorre antes da diagonal principal
      y[i] -= L[i][j] * y[j];       //realiza subtrações de Lx de cada posição após a diagonal principal y = b[i] - L[i][i+1]*y[i] - ... - L[i][n]*y[i]
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
void CalculaXFROMUY(real_t **U, real_t* y, int n, real_t *x){
  for (int i = n-1; i >= 0; i--) {    //percorre diagonal principal de baixo pra cima
    x[i] = y[i];                      //inicia o x com o y da linha
    for (int j = i+1; j < n; j++)     //j inicia após a diagonal principal 
      x[i] -= U[i][j] * x[j];         //realiza subtrações de Uy de cada posição após a diagonal principal x = y[i] - U[i][i+1]*x[i] - ... - U[i][n]*x[i]
    x[i] /= U[i][i];                  //divide pelo elemento da diagonal principal de U[i][i]
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
  \param tTotalY médio para cáculo de y.
  \param tTotalX médio para cáculo de x.
*/
int calculaInversa(real_t **L, real_t **U, real_t **I, real_t *x, real_t *y, int *LUT, unsigned int n, double *tTotalY, double *tTotalX) {

  double tempo;
  int erro;
  real_t *b = (real_t*)malloc(n * sizeof(real_t));
  for(int i = 0; i < n; i++) {
    //inicializa b
    for(int j = 0; j < n; j++) {    
      b[j] = 0.0;
    }
    b[LUT[i]] = 1.0;

    tempo = timestamp();
    CalculaYFROML(L, n, LUT, i, y,b);
    *tTotalY += timestamp() - tempo;

    tempo = timestamp();
    CalculaXFROMUY(U, y, n, x);
    *tTotalX += timestamp() - tempo;

    for(int j = 0; j < n; j++) {
      I[LUT[i]][j] = x[j];
    }
  }

  free(b);
  *tTotalY /= n;
  *tTotalX /= n;
  return 0;
}

int calculaW(real_t **L, real_t **U, real_t **R, real_t **W, real_t *vet, real_t *x, real_t *y, int *LUT, unsigned int n) {

  for(int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++){
      vet[LUT[j]] = R[i][j];
    }

    CalculaYFROML(L, n, LUT, i, y, vet);
    CalculaXFROMUY(U, y, n, x);

    for(int j = 0; j < n; j++) {
      W[LUT[i]][j] = x[j];
    }
  }
  return 0;
}