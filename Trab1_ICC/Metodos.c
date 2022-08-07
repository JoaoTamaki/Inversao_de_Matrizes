#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "utils.h"
#include "sislin.h"
#include "Metodos.h"

//----------------------------------------FUNCOES LAB1----------------------------------------// 

int encontraMax(real_t **M, int *LUT, int i, int n){
  
  int maior = M[LUT[i]][i];
  int indice = LUT[i];
  for (int k = i+1; k < n; k++){
    if (maior < M[LUT[k]][i]){
      maior = M[LUT[k]][i];
      indice = LUT[k];
    }
  }
  return indice;
}

/*void copiaVetor(real_t *a, real_t *b, int k){

  for(int i = 0; i < k; i++){
    b[i] = a[i];
  }
}*/

void retrossubs(SistLinear_t *SL, int *LUT, real_t **x, int n) {  //resolução de SL triangulares

  for (int i = SL->n-1; i >= 0; --i) {
    x[LUT[i]][n] = (SL->b[LUT[i]][i]);                   
    for (int j = i+1; j < SL->n; ++j)
      x[LUT[i]][n] -= SL->A[LUT[i]][j] * x[LUT[j]][n];
    x[LUT[i]][n] /= SL->A[i][i];
  }
}

/*!
  \brief Método da Fatoração LU. Decompõe a Matriz A em L(na parte inferior da diagonal principal) e U(da diagonal
         principal para cima)

  \param SL Ponteiro para o sistema linear
  \param LUT ponteiro para a look up table
  \param tTotal tempo total em milisegundos gastos pelo método
  \param n numero da linha de b
  
  \return código de erro. 0 em caso de sucesso.
*/
int fatoraLU (SistLinear_t *SL, int *LUT, double *tTotal)
{
  //*tTotal = timestamp();
  /* para cada linha a partir da primeira */
  
  for (int i=0; i < SL->n; ++i) {
    int iPivo = encontraMax(SL->A, LUT, i, SL->n);
    if ( LUT[i] != iPivo ){
      int aux = LUT[i];
      LUT[i] = iPivo;
      LUT[iPivo] = aux;
    }

    for(int k=i+1; k < SL->n; ++k) {
      double m = SL->A[LUT[k]][i] / SL->A[LUT[i]][i];           //Calcula multiplicador
      SL->A[LUT[k]][i] = m;                                     //Atualiza L
      for(int j=i+1; j < SL->n; ++j){
        SL->A[LUT[k]][j] -= SL->A[LUT[i]][j] * m;               //Atualiza U
      }
      for(int j=0; j < SL->n; j++){                             //Atualiza a matriz identidade
        if (SL->b[LUT[i]][j] != 0.0)                            //Só faz quando a posição i,j não for 0 
          SL->b[LUT[k]][j] -= SL->b[LUT[i]][j] * m;             //calculo do novo b
      }      
    }
  }

  //*tTotal = timestamp() - *tTotal;
  return 0;
}

/*!
  \brief Método da Eliminação de Gauss

  \param SL Ponteiro para o sistema linear
  \param x ponteiro para o vetor solução
  \param tTotal tempo total em milisegundos gastos pelo método
  \param n numero da linha de b

  \return código de erro. 0 em caso de sucesso.
*/
//int eliminacaoGauss (SistLinear_t *SL, int *LUT, real_t *x, double *tTotal)
//{
  //*tTotal = timestamp();
  /* para cada linha a partir da primeira */
/*  
  for (int i=0; i < SL->n; ++i) {
    int iPivo = encontraMax(SL->A, LUT, i, SL->n);
    //printf("%d\n", iPivo);
    if ( i != iPivo ){
      int aux = LUT[i];
      LUT[i] = iPivo;
      LUT[iPivo] = aux;
    }
    for(int k=i+1; k < SL->n; ++k) {
      double m = SL->A[k][i] / SL->A[i][i];
      SL->A[k][i] = 0.0;
      for(int j=i+1; j < SL->n; ++j)
        SL->A[k][j] -= SL->A[i][j] * m;
      for(int j=0; j < SL->n; j++){                       //Atualiza a matriz identidade
        if (SL->b[i][j] != 0.0)                           //Só faz quando a posição i,j não for 0 
          SL->b[k][j] -= SL->b[i][j] * m;                 //Calculo do novo b
      }  
    }
  }
  //retrossubs(SL, x);
  //*tTotal = timestamp() - *tTotal;
  return 0;
}*/

/*!
  \brief Essa função calcula a norma L2 do resíduo de um sistema linear 

  \param SL Ponteiro para o sistema linear
  \param x Solução do sistema linear
*/
/*real_t normaL2Residuo(SistLinear_t *SL, real_t *x)
{
  real_t *residuo = (real_t *) malloc(SL->n * sizeof(real_t));
  real_t *Ax = (real_t *) malloc(SL->n * sizeof(real_t));
  
  for (int i = 0; i < SL->n; i++){
    for (int j = 0; j < SL->n; j++){
      Ax[i] = SL->A[i][j] * x[j];
    }
    residuo[i] = SL->b[i] - Ax[i];
  }
  free(Ax);
    
  real_t sum = 0.0;
  for(int i = 0; i < SL->n - 1; i++){
    sum += residuo[i] * residuo[i];
  }
    
  free(residuo);

  return sqrt(sum);
}*/

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

//----------------------------------------FUNCOES AUX----------------------------------------// 

int encontraMaxColunaPivo(double** M, int pivNum, int n){
  int maxIndx = pivNum;
  double max = fabs(M[pivNum][pivNum]);
  for(int i= pivNum+1; i<n; i++){
    if(fabs(M[i][pivNum])> max){
      max = fabs(M[i][pivNum]);
      maxIndx = i;
    }
  }
  return maxIndx;
}

void criaMatrizIdentidade(real_t **M, int n){
    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            if(i==j){
                M[i][j]=1.0;
            }
            else{
                M[i][j]=0.0;
            }
        }
    }
    return;
}

/*!
  \brief Essa função multiplicas matrizes quadradas de mesma ordem. Ou seja, mR = mA x mB

  \param mA matriz da esquerda participante da multiplicação. 
  \param mB matriz da direita participante da multiplicação.
  \param mR matriz da resultante da multiplicação.
  \param n  ordem das matrizes quadradas.

  \return 

*/
void MultiplicaMQs(double** mA, double** mB, double** mR, int n){
    for(int mri=0; mri<n; mri++){
        for(int mrj=0; mrj<n; mrj++){
            mR[mri][mrj] = 0.0;
            for(int i=0; i<n; i++){
                mR[mri][mrj] += mA[mri][i] * mB[i][mrj];
            }
        }
    }
}

//----------------------------------------FUNCOES LU----------------------------------------// 
/*!
  \brief Fatoração LU. Decompõe a matriz U em L, U e P que guarda a ordem das trocas de linha. 

  \param L matriz triangular inferior calculada com diagonal principal igual a 1.
  \param U matriz triangular superior calculada.
  \param P vetor que guarda as trocas de linhas durante o pivoteamento. 
  \param n ordem das matrizes quadradas.

  \return êxito da função

*/
int FatoracaoLUMQCOMPIVO(double** L, double** U, int n, int* P){//Talvez void não sei o q pode dar errado nessa funcao

    for(int piv=0; piv<n; piv++){
        
        int newPiv = encontraMaxColunaPivo(U, piv, n);

        if(piv!=newPiv){
            //Troca Linhas 
            trocaLinhasMQ(U, n, piv, newPiv);
            trocaLinhasMQ(L, n, piv, newPiv);
            TrocaElementosVetor(P, piv, newPiv);
        }        

        L[piv][piv] = 1.0;
        double mPiv = U[piv][piv];
        for(int i=piv+1; i<n;i++){
            float mi = U[i][piv];
            U[i][piv] = 0.0;
            L[i][piv] = (mi/mPiv); //não se espera divisão por zero devido pivoteamento
            for(int j=piv+1; j<n; j++){
                U[i][j] = U[i][j] - (mi/mPiv) *U[piv][j];
            }
        }
    }

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


//AINDA NÃO ESTAMOS USANDO -> SERÁ UTILIZADA NA MAIN NA PARTE QUE ESTÁ COMENTADA -> TEM QUE MUDAR (SL_COPIA->A) PARA L E U

/*!
  \brief Cácula Y da equação "L Y = I" como parte da fatoração LU.
         Utiliza várias (forward) substituições onde uma coluna de x "Yci" pode ser calculada através de L e uma coluna de "I";
         Ou seja, ao invés de fazer a (forward) substituição da forma convencional Ay=b. É feita coluna a coluna até termos todas colunas de Y.
         Para simplificar essa conta Y é considerado já iniciado como I.
         Além disso, as divisões que apareceriam não aparecem, pois 1 é o elemento neutro da multiplicação.
         E, a matriz L contem somente somente "1"s na sua diagonal principal.  

  \param Y matriz a ser cáculada previamente iniciada como I.
  \param L matriz da fatoração LU. Matriz triangular inferior cuja diagonal principal consiste em somente 1's. 

*/
void CalculaYFROML(SistLinear_t *SL, int *LUT, double** Y){
    for(int yj = 0; yj < SL->n; yj++){
        for(int i = 0; i < SL->n; i++){
            for(int j = 0; j < i; j++){
                Y[LUT[i]][yj] -= SL->A[LUT[i]][j] * Y[LUT[j]][yj];
            }
        }
    }
    return;
}

/*!
  \brief Cácula X da equação "U X = Y" como parte da fatoração LU.
         Utiliza várias retrosubstiruições onde uma coluna de x "Xci" pode ser calculada através de U e uma coluna de "Yci";
         Ou seja, ao invés de fazer a retrosubstituição da forma convencional Ax=b. É feita coluna a coluna até termos todas colunas de X.

  \param X matriz a ser cáculada.
  \param U matriz da fatoração LU. Matriz triangular superior. 
  \param Y matriz da previamente calculada usando L.
  \param n  ordem das matrizes quadradas.

*/
void CalculaXFROMUY(double** X, double** U, double** Y, int n){
    for(int xj=0; xj<n; xj++){  
        //Retrosubstituicao Adaptada
        for(int i= n-1; i>=0; i--){
            X[i][xj]=Y[i][xj];
            for(int j= i+1; j<n; j++){
                X[i][xj] -= U[i][j] * X[j][xj];
            }
            X[i][xj] /= U[i][i];
        }
    }
}

real_t calculaDeterminante(real_t **U, int n) {

  real_t determinante = 1.0;

  for (int i = 0;  i < n; ++i) 
    determinante = determinante * U[i][i];
  return determinante;
}