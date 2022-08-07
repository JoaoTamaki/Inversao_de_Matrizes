#ifndef __METODOS_H__
#define __METODOS_H__

// Parâmetros para teste de convergência
#define MAXIT 100     // número máximo de iterações para métodos iterativos
#define ERRO 1.0e-6   // Tolerância para critérios de parada em métodos iterativos

//----------------------------------------FUNCOES LAB1----------------------------------------// 
int encontraMax(real_t **M, int *LUT, int i, int n);
//void copiaVetor(real_t *a, real_t *b, int k);

void retrossubs(SistLinear_t *SL, int *LUT, real_t **x, int n);
//int eliminacaoGauss (SistLinear_t *SL, int *LUT, real_t *x, double *tTotal);

//real_t normaL2Residuo(SistLinear_t *SL, real_t *x);
//int refinamento (SistLinear_t *SL, real_t *x, real_t erro, double *tTotal, double *norma_ref);

//----------------------------------------FUNCOES AUX----------------------------------------// 
int encontraMaxColunaPivo(double** M, int pivNum, int n);
void criaMatrizIdentidade(real_t **M, int n);
void MultiplicaMQs(double** mA, double** mB, double** mR, int n);

//----------------------------------------FUNCOES LU----------------------------------------// 
int FatoracaoLUMQCOMPIVO(double** L, double** U, int n, int* P);

void TrocaElementosVetor(int *vet, int i1, int i2);
void trocaLinhasMQ(double **M, int n, int l1, int l2);

int fatoraLU (SistLinear_t *SL, int *LUT, double *tTotal);
void CalculaYFROML(SistLinear_t *SL, int *LUT, double** Y);
void CalculaXFROMUY(double** X, double** U, double** Y, int n);
real_t calculaDeterminante(real_t **U, int n);

#endif // __METODOS_H__