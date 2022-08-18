#ifndef __METODOS_H__
#define __METODOS_H__

// Parâmetros para teste de convergência
#define MAXIT 100     // número máximo de iterações para métodos iterativos
#define ERRO 1.0e-6   // Tolerância para critérios de parada em métodos iterativos

//----------------------------------------FUNCOES AUX----------------------------------------// 
void calculaNovoI(real_t **I, real_t **W, unsigned int n);
real_t normaL2Residuo(real_t **R, unsigned int n);
real_t refinamento(SistLinear_t *SL, real_t **L, real_t **U, real_t **I, real_t **R, real_t **W, real_t *vet, real_t *x, real_t *y, int *LUT, double *tTempoIter, double *tTempoResiduo);
void calculaMatrizResiduo(real_t** mA, real_t** mB, real_t** mI, real_t** mR, int n);
int encontraMaxColunaPivo(double** M, int pivNum, int n);
void criaMatrizIdentidade(real_t **M, int n);
real_t calculaDeterminante(real_t **M, int n);

//----------------------------------------FUNCOES LU----------------------------------------// 
int FatoracaoLU_PivoParcial(double** L, double** U, int n, int* P, double *tTotal);
void TrocaElementosVetor(int *vet, int i1, int i2);
void trocaLinhasMQ(double **M, int n, int l1, int l2);
void CalculaYFROML(real_t **L, int n, int *LUT, int k, real_t* y, real_t* b);
void CalculaXFROMUY(real_t **U, real_t* y, int n, real_t *x);
int calculaInversa(real_t **L, real_t **U, real_t **I, real_t *x, real_t *y, int *LUT, unsigned int n, double *tTotalY, double *tTotalX);
int calculaW(real_t **L, real_t **U, real_t **R, real_t **W, real_t *vet, real_t *x, real_t *y, int *LUT, unsigned int n);

#endif // __METODOS_H__