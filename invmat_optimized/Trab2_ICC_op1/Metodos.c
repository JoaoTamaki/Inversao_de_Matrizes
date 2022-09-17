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
  \param pad Tamanho do pad
  \return
*/
void calculaNovoI(real_t* restrict I, real_t* restrict W, unsigned int n, unsigned int pad){
  int i, j;
  unsigned int mult, np = n + pad;
  for (i = 0; i < n; i++){
    mult = i*np;
    for (j = 0; j < n; j++){
      I[mult+j] += W[mult+j]; 
    }
  }
}

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 
  \param R    Ponteiro para a matriz de Resíduo
  \param n    Tamanho da dimensão da matriz
  \param pad  Tamanho do pad
  \return Norma do resíduo
*/
real_t normaL2Residuo(real_t* R, unsigned int n, unsigned int pad){
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
  \param pad            Tamanho do pad
  \param tTempoIter     Ponteiro para o endereço do tempo de cada iteração (op1): LUx = b -> Ly = b, Ux = y
  \param tTempoResiduo  Ponteiro para o endereço do tempo de cada calculo da matriz Resíduo (op2): R = b - A * A^(-1)
  \return Norma do resíduo
*/
real_t refinamento(SistLinear_t* SL, real_t* L, real_t* U, real_t* I, real_t* R, real_t* W, real_t* vet, real_t* x, real_t* y, int* LUT, unsigned int pad, double* tTempoIter, double* tTempoResiduo)
{
  //Inicializa variáveis
  unsigned int n = SL->n;
  real_t normaResiduo;
  //Realiza Resíduo = B - A * I -> B(Identidade), A(Matriz de entrada), I(Matriz Inversa)
  *tTempoResiduo = timestamp();
  calculaMatrizResiduo(SL->A, SL->b, I, R, n, pad);
  *tTempoResiduo = timestamp() - *tTempoResiduo;
  //Aw = R
  *tTempoIter = timestamp();  
  calculaW(L, U, R, W, vet, x, y, LUT, n, pad);
  *tTempoIter = timestamp() - *tTempoIter;
  //Calcula nova Inversa
  calculaNovoI(I, W, n, pad);
  //Calcula norma do resíduo
  normaResiduo = normaL2Residuo(R, n, pad);
  return normaResiduo;
}

/*!
  \brief Esta função realiza o cálculo da matriz Resíduo (op2): R = b - A * A^(-1)
  \param mA   Ponteiro para a matriz A participante da multiplicação
  \param mB   Ponteiro para a matriz B participante da subtração
  \param mI   Ponteiro para a matriz Inversa participante da multiplicação
  \param mR   Ponteiro para a matriz Resíduo resultante
  \param n    Tamanho da dimensão da matriz
  \param pad  Tamanho do pad
  \return 
*/
void inline calculaMatrizResiduo(real_t* restrict mA, real_t* restrict mB, real_t* restrict mI, real_t* restrict mR, unsigned int n, unsigned int pad) {
  //LIKWID_MARKER_INIT;
  //LIKWID_MARKER_START("op2");
  
  unsigned int multj, mult1, mult2, mult3, mult, np = n + pad;
  int aux1, aux2, aux3;

  for (int i=0; i < n - n % TAMROLL; i += TAMROLL){
    mult = i*np;
    aux1 = i+1;
    mult1 = aux1*np;
    aux2 = i+2;
    mult2 = aux2*np;
    aux3 = i+3;
    mult3 = aux3*np;
    for (int j = 0; j < n; j++){
      mR[mult+j] = 0.0;
      mR[mult1+j] = 0.0;
      mR[mult2+j] = 0.0;
      mR[mult3+j] = 0.0; 
      multj = j*np;
      for (int k = 0; k < n; k++){
        mR[mult+j] += mA[mult+k] * mI[multj+k];
        mR[mult1+j] += mA[mult1+k] * mI[multj+k];
        mR[mult2+j] += mA[mult2+k] * mI[multj+k];
        mR[mult3+j] += mA[mult3+k] * mI[multj+k];
      }
      mR[mult+j] = mB[mult+j] - mR[mult+j];
      mR[mult1+j] = mB[mult1+j] - mR[mult1+j]; 
      mR[mult2+j] = mB[mult2+j] - mR[mult2+j]; 
      mR[mult3+j] = mB[mult3+j] - mR[mult3+j]; 
    }
  }
  for (int i = n - n % TAMROLL; i < n; i++){
    mult = i*np;               //em vez de executar i*j vezes, faz apenas i vezes
    for (int j = 0; j < n; j++) {
      mR[mult+j] = 0.0;
      multj = j*np;
      for (int k = 0; k < n; k++)
        mR[mult+j] += mA[mult+k] * mI[multj+k];
    }
    mR[mult+i] = mB[mult+i] - mR[mult+i]; 
  }
  
  /*!Versão sem unloop and jam 
  unsigned int mult, np = n + pad;
  for (int i = 0; i < n; i ++) {
    mult = i*np;               //em vez de executar i*j vezes, faz apenas i vezes
    for (int j = 0; j < n; j++) {
      mR[mult+j] = 0.0;
      for (int k = 0; k < n; k++)
        mR[mult+j] += mA[mult+k] * mI[j*np+k];
    }
    mR[mult+i] = mB[mult+i] - mR[mult+i]; 
  }
  */
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
unsigned int encontraMaxColunaPivo(real_t* M, unsigned int pivNum, unsigned int n, unsigned int pad) {

  int i;
  unsigned int np = n + pad;
  unsigned int maxIndx = pivNum;
  real_t max = fabs(M[pivNum*np+pivNum]);
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
  \param M Ponteiro para a matriz
  \param n Tamanho da dimensão da matriz
  \param pad Tamanho do pad
  \return
  */
void criaMatrizIdentidade(real_t* M, unsigned int n, unsigned int pad) {

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
  \param M    Ponteiro para a matriz
  \param n    Tamanho da dimensão da matriz
  \param pad  Tamanho do pad
  \return Valor do determinante
*/
real_t calculaDeterminante(real_t* M, unsigned int n, unsigned int pad) {

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
  \param L      Matriz triangular inferior calculada com diagonal principal igual a 1
  \param U      Matriz triangular superior calculada
  \param LUT    Ponteiro para o vetor da Look up table, que contem a ordem das linhas das matrizes L e U do pivoteamento parcial
  \param n      Tamanho da dimensão da matriz
  \param pad    Tamanho do pad
  \param tTotal Tempo total em milisegundos gastos pelo método
  \return Código de erro. 0 em caso de sucesso.
*/
int FatoracaoLU_PivoParcial(real_t* restrict L, real_t* restrict U, unsigned int n, unsigned int pad, int* restrict LUT, double* restrict tTotal) {

  int i, j;
  unsigned int piv, newPiv, np = n + pad;
  real_t mPiv, mi;
  *tTotal = timestamp();
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
  \param pad  Tamanho do pad
  \param l1   Índice do linha 1
  \param l2   Índice do linha 2
  \return 
*/
void trocaLinhasMQ(real_t* M, unsigned int n, unsigned int pad, int l1, int l2){
  real_t aux;
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
  \brief Calcula Y da equação "Ly = b" como parte da fatoração LU.
  \param L    Ponteiro para a matriz L da fatoração LU. Matriz triangular inferior cuja diagonal principal consiste em somente 1's. 
  \param y    Pontero para o vetor y a ser calculado, que é usado para calcula Ux = y.
  \param b    Ponteiro para o vetor b
  \param n    Tamanho da dimensão da matriz
  \param pad  Tamanho do pad
  \return 
*/
void CalculaYFROML(real_t* restrict L, real_t* restrict y, real_t* restrict b, unsigned int n, unsigned int pad) {
  unsigned int mult, mult1, mult2, mult3, aux1, aux2, aux3, np = n + pad;
  
  for (int i=0; i < n - n % TAMROLL; i += TAMROLL) {    //percorre diagonal principal de cima pra baixo
    //inicializações -> inicia o y com o b da linha, calcula mult para em vez de executar i*j vezes, faz apenas i vezes e calcula auxiliares para auxiliar no calculo da linha e economizar operações 
    mult = i*np;
    y[i] = b[i];
    aux1 = i+1;
    mult1 = aux1*np;
    y[aux1] = b[aux1];
    aux2 = i+2;
    mult2 = aux2*np;
    y[aux2] = b[aux2];
    aux3 = i+3;
    mult3 = aux3*np;
    y[aux3] = b[aux3];
    //realiza conta do bloco de [i,0] até o [i+3,i]
    for (int j = 0; j < i; j++){        //loop fusion da parte do bloco  -> //realiza subtrações de Lx de cada posição após a diagonal principal y = b[i] - L[i][i+1]*y[i] - ... - L[i][n]*y[i]    
      y[i] -= L[mult+j] * y[j];
      y[aux1] -= L[mult1+j] * y[j]; 
      y[aux2] -= L[mult2+j] * y[j]; 
      y[aux3] -= L[mult3+j] * y[j]; 
    }
    //realiza conta da parte inferior escada do [i,i] até [i+3,i+3] (menos na diagonal principal, pois tudo dividido por 1 da ele mesmo)
    y[aux1] -= L[mult1+i] * y[i];
    y[aux2] -= L[mult2+i] * y[i];
    y[aux3] -= L[mult3+i] * y[i];
    y[aux2] -= L[mult2+aux1] * y[aux1];
    y[aux3] -= L[mult3+aux1] * y[aux1];
    y[aux3] -= L[mult3+aux2] * y[aux2];           
  }
  for (int i = n - n % TAMROLL; i < n; i++){            //percorre diagonal principal de cima pra baixo
    y[i] = b[i];                                        //inicia o y com o b da linha
    mult = i*np;                                        //em vez de executar i*j vezes, faz apenas i vezes
    for (int j = 0; j < i; j++)                             //j inicia no inicio e percorre antes da diagonal principal
      y[i] -= L[mult+j] * y[j];                         //realiza subtrações de Lx de cada posição após a diagonal principal y = b[i] - L[i][i+1]*y[i] - ... - L[i][n]*y[i]
  }

  /*! Versão sem loop unroll & jam
  int i, j;
  unsigned int mult, np = n + pad;
  for (i = 0; i < n; i++) {               //percorre diagonal principal de cima pra baixo
    y[i] = b[i];                          //inicia o y com o b da linha
    mult = i*np;                          //em vez de executar i*j vezes, faz apenas i vezes
    for (j = 0; j < i; j++)               //j inicia no inicio e percorre antes da diagonal principal
      y[i] -= L[mult+j] * y[j];           //realiza subtrações de Lx de cada posição após a diagonal principal y = b[i] - L[i][i+1]*y[i] - ... - L[i][n]*y[i]
  }*/
  return;
}

/*!
  \brief Calcula X da equação "Ux = Y" como parte da fatoração LU.
  \param U    Ponteiro para a matriz da fatoração LU. Matriz triangular superior.
  \param x    Ponteiro para o vetor x a ser calculado, que é uma linha da matriz identidade.
  \param y    Ponteiro para o vetor y previamente calculada na função CalculaYFROML.
  \param n    Tamanho da dimensão da matriz
  \param pad  Tamanho do pad
*/
void CalculaXFROMUY(real_t* restrict U, real_t* restrict x, real_t* restrict y, unsigned int n, unsigned int pad){

  unsigned int mult, mult1, mult2, mult3, aux1, aux2, aux3, np = n + pad;

  for (int i = n-1; i >= TAMROLL - 1; i -= TAMROLL) {    //percorre diagonal principal de cima pra baixo
    //inicializações -> inicia o y com o b da linha, calcula mult para em vez de executar i*j vezes, faz apenas i vezes e calcula auxiliares para auxiliar no calculo da linha e economizar operações 
    aux3 = i-3;
    mult3 = aux3*np;
    x[aux3] = y[aux3];
    aux2 = i-2;
    mult2 = aux2*np;
    x[aux2] = y[aux2];
    aux1 = i-1;
    mult1 = aux1*np;
    x[aux1] = y[aux1];
    mult = i*np;
    x[i] = y[i];

    //realiza conta do bloco de [i-3,0] até o [i,i]
    for (int j = i+1; j < n; j++){        //loop fusion da parte do bloco  -> realiza subtrações de Uy de cada posição após a diagonal principal x = y[i] - U[i][i+1]*x[i] - ... - U[i][n]*x[i]
      x[aux3] -= U[mult3+j] * x[j]; 
      x[aux2] -= U[mult2+j] * x[j];
      x[aux1] -= U[mult1+j] * x[j];
      x[i] -= U[mult+j] * x[j];
    }        
    //realiza conta da parte inferior escada do [i,i-3] até [i+3,i] (menos na diagonal principal, que será dividida)

    x[i] /= U[mult+i];                            //divide pelo elemento da diagonal principal de U[i][i]

    x[aux3] -= U[mult3+i] * x[i];
    x[aux2] -= U[mult2+i] * x[i];
    x[aux1] -= U[mult1+i] * x[i];
    x[aux1] /= U[mult1+aux1];                    //divide pelo elemento da diagonal principal de U[i][i]
    
    x[aux3] -= U[mult3+aux1] * x[aux1];
    x[aux2] -= U[mult2+aux1] * x[aux1];
    x[aux2] /= U[mult2+aux2];                    //divide pelo elemento da diagonal principal de U[i][i]
    
    x[aux3] -= U[mult3+aux2] * x[aux2];
    x[aux3] /= U[mult3+aux3];                    //divide pelo elemento da diagonal principal de U[i][i]

  }

  for (int i = (n % TAMROLL) - 1; i >= 0; i--){   //percorre diagonal principal de baixo pra cima
    x[i] = y[i];                                  //inicia o x com o y da linha
    mult = i*np;                                  //em vez de executar i*(np-j+1)+1 vezes, faz apenas i+1 vezes
    for (int j = i+1; j < n; ++j){                //j inicia após a diagonal principal 
      x[i] -= U[mult+j] * x[j];                   //realiza subtrações de Uy de cada posição após a diagonal principal x = y[i] - U[i][i+1]*x[i] - ... - U[i][n]*x[i]
    }
    x[i] /= U[mult+i];                            //divide pelo elemento da diagonal principal de U[i][i]
  }


  /*! Versão sem loop unroll & jam
  int i, j;
  unsigned int mult, np = n + pad;
  for (i = n-1; i >= 0; --i) {            //percorre diagonal principal de baixo pra cima
    x[i] = y[i];                          //inicia o x com o y da linha
    mult = i*np;                          //em vez de executar i*(np-j+1)+1 vezes, faz apenas i+1 vezes
    for (j = i+1; j < n; ++j){            //j inicia após a diagonal principal 
      x[i] -= U[mult+j] * x[j];           //realiza subtrações de Uy de cada posição após a diagonal principal x = y[i] - U[i][i+1]*x[i] - ... - U[i][n]*x[i]
    }
    x[i] /= U[mult+i];                    //divide pelo elemento da diagonal principal de U[i][i]
  }*/
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
  \param pad  Tamanho do pad
  \return
*/
void calculaInversa(real_t* L, real_t* U, real_t* I, real_t* x, real_t* y, int* LUT, unsigned int n, unsigned int pad) {
  unsigned int np = n + pad;
  real_t *b = aligned_alloc(16, (np)*sizeof(real_t));
  //LIKWID_MARKER_INIT; 
  //LIKWID_MARKER_START("op1");
  for(int i = 0; i < n; i++) {
    //inicializa b
    for(int j = 0; j < n; j++) {    
      b[j] = 0.0;
    }
    b[LUT[i]] = 1.0;
    CalculaYFROML(L, y, b, n, pad);
    CalculaXFROMUY(U, x, y, n, pad);
    for(int j = 0; j < n; j++) {
      I[LUT[i]*np+j] = x[j];
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
  \param pad  Tamanho do pad.
  \return
*/
void calculaW(real_t* L, real_t* U, real_t* R, real_t* W, real_t* vet, real_t* x, real_t* y, int* LUT, unsigned int n, unsigned int pad) {
  int i, j;
  unsigned int mult, np = n + pad;
  LIKWID_MARKER_INIT;
  LIKWID_MARKER_START("op1");
  for(i = 0; i < n; i++) {
    mult = i*np;                          //em vez de executar i*j vezes, faz apenas i vezes
    //Invete R para usar no LUw=R
    for (j = 0; j < n; j++){
      vet[LUT[j]] = R[mult+j];
    }
    CalculaYFROML(L, y, vet, n, pad);
    CalculaXFROMUY(U, x, y, n, pad);
    mult = LUT[i]*np;
    for(j = 0; j < n; j++) {
      W[mult+j] = x[j];
    }
  }
  LIKWID_MARKER_STOP("op1");
  LIKWID_MARKER_CLOSE;
}