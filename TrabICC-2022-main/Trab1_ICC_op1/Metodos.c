#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

//----------------------------------------FUNÇÕES AUX----------------------------------------// 
/*!
  \brief Essa função calcula a nova inversa ataulizada pelo refinamento
  \param I Ponteiro para a matriz Inversa
  \param W Ponteiro para a matriz W
  \param n Tamanho da dimensão da matriz
  \return
*/
void calculaNovoI(real_t** I, real_t** W, unsigned int n){
  
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      I[i][j] += W[i][j]; 
    }
  }
}

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 
  \param R Ponteiro para a matriz de resíduo
  \param n Tamanho da dimensão da matriz
  \return Norma do resíduo
*/
real_t normaL2Residuo(real_t** R, unsigned int n){
  
  real_t sum;
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n; j++){
      sum += R[i][j]*R[i][j];
    }
  }
  return sqrt(sum);
}

/*!
  \brief Essa função realiza o refinamento
  \param SL             Ponteiro para o sistema linear
  \param L              Ponteiro para a matriz L da Fatoração LU
  \param U              Ponteiro para a matriz U da Fatoração LU
  \param I              Ponteiro para a matriz Inversa
  \param R              Ponteiro para a matriz de Resíduo
  \param W              Ponteiro para a matriz W
  \param vet            Ponteiro para o vetor auxiliar vet, utilizado no calculo de W para inverter as linhas do resíduo
  \param x              Ponteiro para o vetor auxiliar x, utilizado no calculo de W para o calculo de x de cada linha da matriz
  \param y              Ponteiro para o vetor auxiliar x, utilizado no calculo de W para o calculo de y de cada linha da matriz
  \param LUT            Ponteiro para o vetor da Look up table, que contem a ordem das linhas das matrizes L e U do pivoteamento parcial
  \param tTempoIter     Ponteiro para o endereço do tempo de cada iteração (op1): LUx = b -> Ly = b, Ux = y
  \param tTempoResiduo  Ponteiro para o endereço do tempo de cada calculo da matriz Resíduo (op2): R = b - A * A^(-1)
  \return Norma do resíduo
*/
real_t refinamento(SistLinear_t* SL, real_t** L, real_t** U, real_t** I, real_t** R, real_t** W, real_t* vet, real_t* x, real_t* y, int* LUT, double* tTempoIter, double* tTempoResiduo){
  
  //Inicializa variáveis
  unsigned int n = SL->n;
  real_t normaResiduo;
  //Realiza Resíduo = B - A * I -> B(Identidade), A(Matriz de entrada), I(Matriz Inversa)
  *tTempoResiduo = timestamp();
  calculaMatrizResiduo(SL->A, SL->b, I, R, n);
  *tTempoResiduo = timestamp() - *tTempoResiduo;
  //Aw = R
  *tTempoIter = timestamp();
  calculaW(L, U, R, W, vet, x, y, LUT, n);
  *tTempoIter = timestamp() - *tTempoIter;
  //Calcula nova Inversa
  calculaNovoI(I, W, n);
  //Calcula norma do resíduo
  normaResiduo = normaL2Residuo(R, n);
  return normaResiduo;
}

/*!
  \brief Esta função realiza o cálculo da matriz Resíduo (op2): R = b - A * A^(-1)
  \param mA   Ponteiro para a matriz A participante da multiplicação
  \param mB   Ponteiro para a matriz B participante da subtração
  \param mI   Ponteiro para a matriz Inversa participante da multiplicação
  \param mR   Ponteiro para a matriz Resíduo resultante
  \param n    Tamanho da dimensão da matriz
  \return 
*/
void calculaMatrizResiduo(real_t** mA, real_t** mB, real_t** mI, real_t** mR, int n) {
  //LIKWID_MARKER_INIT;
  //LIKWID_MARKER_START("op2");
  for (int i = 0; i < n; i ++) {
    for (int j = 0; j < n; j++) {
      mR[i][j] = 0.0;
      for (int k = 0; k < n; k++)
        mR[i][j] += mA[i][k] * mI[j][k];
    }
    mR[i][i] = mB[i][i] - mR[i][i]; 
  }
  //LIKWID_MARKER_STOP("op2");
  //LIKWID_MARKER_CLOSE;
  return;  
}

/*!
  \brief Encontra o maior valor de uma coluna: pivô
  \param M      Ponteiro para a matriz
  \param pivNum Índice da linha em análise que pode ser trocada pela linha pivô
  \param n      Tamanho da dimensão da matriz
  \param pad    Tamanho do pad
  \return Índice da linha correspondente ao pivô (que tem o maior valor da coluna)
*/
int encontraMaxColunaPivo(real_t** M, int pivNum, int n) {

  int maxIndx = pivNum;
  real_t max = fabs(M[pivNum][pivNum]);

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
  \param M Ponteiro para a matriz
  \param n Tamanho da matriz
  \return 
  */
void criaMatrizIdentidade(real_t** M, int n) {

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
  \param M Ponteiro para a matriz
  \param n Tamanho da dimensão da matriz
  \return Valor do determinante
*/
real_t calculaDeterminante(real_t** M, int n) {

  real_t determinante = 1.0;
  for (int i = 0;  i < n; ++i) 
    determinante = determinante * M[i][i];
  return determinante;
}

//----------------------------------------FUNÇÕES LU----------------------------------------// 
/*!
  \brief Fatoração LU. Decompõe a matriz U em L, U e P que guarda a ordem das trocas de linha. Para isso, a função executa uma Eliminação de Gauss com pivoteamento parcial. 
  \param L      Matriz triangular inferior calculada com diagonal principal igual a 1
  \param U      Matriz triangular superior calculada
  \param LUT    Ponteiro para o vetor da Look up table, que contem a ordem das linhas das matrizes L e U do pivoteamento parcial 
  \param n      Tamanho da dimensão da matriz
  \param tTotal Tempo total em milisegundos gastos pelo método
  \return Código de erro. 0 em caso de sucesso.
*/
int FatoracaoLU_PivoParcial(real_t** L, real_t** U, int n, int* LUT, double *tTotal) {

  *tTotal = timestamp();
  // ESCALONAMENTO: PIVOTEAMENTO PARCIAL DA MATRIZ A COPIADA EM U, 
  // COM TRIANGULAÇÃO RESULTANDO EM MATRIZ SUPERIOR
  for (int piv = 0; piv < n; piv++) {               // para cada equação (linha) da matriz A 
    int newPiv = encontraMaxColunaPivo(U, piv, n);  // encontra o maior valor da coluna correspondente
    
    if (piv != newPiv) {                            // se o índice da linha pivô for diferente do índice da equação atual
      trocaLinhasMQ(U, n, piv, newPiv);             // troca linha de U atual pela linha pivô
      trocaLinhasMQ(L, n, piv, newPiv);             // troca linha de L atual pela linha pivô
      TrocaElementosVetor(LUT, piv, newPiv);          // armazena a troca na LUT
    }  
    L[piv][piv] = 1.0;                              // insere 1 na linha atual, na posição da diagonal principal de L
    real_t mPiv = U[piv][piv];

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

/*!
  \brief Função para troca de elementos em um vetor
  \param vet  Ponteiro para o vetor
  \param i1   Índice do elemento 1
  \param i2   Índice do elemento 2
  \return 
*/
void TrocaElementosVetor(int* vet, int i1, int i2){
  int aux = vet[i1];
  vet[i1] = vet[i2];
  vet[i2] = aux;
  return;
}

/*!
  \brief Função para troca de elementos entre vetores
  \param vet  Ponteiro para a matriz
  \param n    Tamanho da dimensão da matriz
  \param l1   Índice do linha 1
  \param l2   Índice do linha 2
  \return 
*/
void trocaLinhasMQ(real_t** M, int n, int l1, int l2){
  real_t* aux;
  aux = M[l1];
  M[l1] = M[l2];
  M[l2] = aux;
  return;
}

/*!
  \brief Calcula Y da equação "Ly = b" como parte da fatoração LU.
  \param L    Ponteiro para a matriz L da fatoração LU. Matriz triangular inferior cuja diagonal principal consiste em somente 1's. 
  \param y    Pontero para o vetor y a ser calculado, que é usado para calcula Ux = y.
  \param b    Ponteiro para o vetor b
  \param n    Tamanho da dimensão da matriz
  \return 
*/
void CalculaYFROML(real_t** L, real_t* y, real_t* b, int n) {
  for (int i = 0; i < n; i++) {     //percorre diagonal principal de cima pra baixo
    y[i] = b[i];                    //inicia o y com o b da linha
    for (int j = 0; j < i; j++)     //j inicia no inicio e percorre antes da diagonal principal
      y[i] -= L[i][j] * y[j];       //realiza subtrações de Lx de cada posição após a diagonal principal y = b[i] - L[i][i+1]*y[i] - ... - L[i][n]*y[i]
  }
  return;
}

/*!
  \brief Calcula X da equação "Ux = Y" como parte da fatoração LU.
  \param U    Ponteiro para a matriz da fatoração LU. Matriz triangular superior.
  \param x    Ponteiro para o vetor x a ser calculado, que é uma linha da matriz identidade.
  \param y    Ponteiro para o vetor y previamente calculada na função CalculaYFROML.
  \param n    Tamanho da dimensão da matriz
*/
void CalculaXFROMUY(real_t** U, real_t* x, real_t* y, int n){
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
  \param L    Ponteiro para a matriz da fatoração LU. Matriz triangular inferior cuja diagonal principal consiste em somente 1's.
  \param U    Ponteiro para a matriz da fatoração LU. Matriz triangular superior.  
  \param I    Ponteiro para a matriz que armazenará a inversa de A.
  \param x    Ponteiro para o vetor x a ser calculado, que é usado para calcula Ux = y.
  \param y    Ponteiro para o vetor y a ser calculado, que é usado para calcula Ly = b.
  \param LUT  Ponteiro para a vetor Look up table que armazena as trocas de linha em b.
  \param n    Tamanho da dimensão da matriz
  \return
*/
void calculaInversa(real_t** L, real_t** U, real_t** I, real_t* x, real_t* y, int* LUT, unsigned int n) {
  //LIKWID_MARKER_INIT;
  //LIKWID_MARKER_START("op1");
  real_t *b = (real_t*)malloc(n * sizeof(real_t));
  for(int i = 0; i < n; i++) {
    //inicializa b
    for(int j = 0; j < n; j++) {    
      b[j] = 0.0;
    }
    b[LUT[i]] = 1.0;
    CalculaYFROML(L, y, b, n);
    CalculaXFROMUY(U, x, y, n);
    for(int j = 0; j < n; j++) {
      I[LUT[i]][j] = x[j];
    }
  }
  //LIKWID_MARKER_STOP("op1");
  //LIKWID_MARKER_CLOSE;
  free(b);
}

/*!
  \brief Calcula matriz inversa usando fatoração LU e eliminação de Gauss.
  \param L    Ponteiro para a matriz da fatoração LU. Matriz triangular inferior cuja diagonal principal consiste em somente 1's.
  \param U    Ponteiro para a matriz da fatoração LU. Matriz triangular superior.  
  \param R    Ponteiro para a matriz que armazenará o Resíduo.
  \param w    Ponteiro para a matriz que armazenará o Resíduo.
  \param vet  Ponteiro para o vetor auxiliar vet, utilizado no calculo de W para inverter as linhas do resíduo.
  \param x    Ponteiro para o vetor x a ser calculado, que é usado para calcula Ux = y.
  \param y    Ponteiro para o vetor y a ser calculado, que é usado para calcula Ly = b.
  \param LUT  Ponteiro para o vetor da Look up table, que contem a ordem das linhas das matrizes L e U do pivoteamento parcial.
  \param n    Tamanho da dimensão da matriz.
  \return
*/
void calculaW(real_t** L, real_t** U, real_t** R, real_t** W, real_t* vet, real_t* x, real_t* y, int* LUT, unsigned int n) {
  LIKWID_MARKER_INIT;
  LIKWID_MARKER_START("op1");
  for(int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++){
      vet[LUT[j]] = R[i][j];
    }
    CalculaYFROML(L, y, vet, n);
    CalculaXFROMUY(U, x, y, n);
    for(int j = 0; j < n; j++) {
      W[LUT[i]][j] = x[j];
    }
  }
  LIKWID_MARKER_STOP("op1");
  LIKWID_MARKER_CLOSE;
}