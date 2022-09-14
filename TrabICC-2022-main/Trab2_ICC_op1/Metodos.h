#ifndef __METODOS_H__
#define __METODOS_H__

// Parâmetros para teste de convergência
#define MAXIT 100     // número máximo de iterações para métodos iterativos
#define ERRO 1.0e-6   // Tolerância para critérios de parada em métodos iterativos
#define TAMROLL 4    // Tamanho do unroll 
//----------------------------------------FUNCOES AUX----------------------------------------// 
void calculaNovoI(real_t* I, real_t* W, unsigned int n, unsigned int pad);
real_t normaL2Residuo(real_t* R, unsigned int n, unsigned int pad);
real_t refinamento(SistLinear_t* SL, real_t* L, real_t* U, real_t* I, real_t* R, real_t* W, real_t* vet, real_t* x, real_t* y, int* LUT, unsigned int pad, double* tTempoIter, double* tTempoResiduo);
void calculaMatrizResiduo(real_t* mA, real_t* mB, real_t* mI, real_t* mR, unsigned int n, unsigned int pad);
unsigned int encontraMaxColunaPivo(real_t* M, unsigned int pivNum, unsigned int n, unsigned int pad);
void criaMatrizIdentidade(real_t* M, unsigned int n, unsigned int pad);
real_t calculaDeterminante(real_t* M, unsigned int n, unsigned int pad);

//----------------------------------------FUNCOES LU----------------------------------------// 
int FatoracaoLU_PivoParcial(real_t* L, real_t* U, unsigned int n, unsigned int pad, int* LUT, double* tTotal);
void TrocaElementosVetor(int* vet, int i1, int i2);
void trocaLinhasMQ(real_t* M, unsigned int n, unsigned int pad, int l1, int l2);
void CalculaYFROML(real_t* L, real_t* y, real_t* b, unsigned int n, unsigned int pad);
void CalculaXFROMUY(real_t* U, real_t* x, real_t* y, unsigned int n, unsigned int pad);
void calculaInversa(real_t* L, real_t* U, real_t* I, real_t* x, real_t* y, int* LUT, unsigned int n, unsigned int pad);
void calculaW(real_t* L, real_t* U, real_t* R, real_t* W, real_t* vet, real_t* x, real_t* y, int* LUT, unsigned int n, unsigned int pad);

#endif // __METODOS_H__