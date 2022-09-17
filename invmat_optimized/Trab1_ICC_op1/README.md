# Trabalho 1 - Introdução à Computação Científica: inversão de matriz usando fatoração LU

## Equipe

- João Victor Kenji Tamaki - GRR20182965 
- Marisa Sel Franco - GRR20186556

## Descrição

O objetivo deste trabalho é implementar um programa computacional que, dada uma matriz quadrada A de dimensão n, **encontre a matriz inversa de A(A−1)**, tal que A×A−1=I, onde I é a matriz identidade. Para tal, o programa usa o Método da Eliminação de Gauss com Pivoteamento Parcial, Fatoração LU e Refinamento:
- Cálculo das matrizes L e U (Eliminação de Gauss com Pivotamento Parcial);
- Solução parcial do Sistema Linear (Retrosubstituição);
- Solução final do Sistema Linear por Refinamento a partir do Resíduo.

A medida de erro a cada iteração do refinamento é a norma L² do resíduo (|| R ||):

    R=I−A∗A−1

||R||= sqrt(∑R[i,j]²),   1 ⩽ i, j ⩽ n

Neste documento, será explicada de forma breve a compilação, execução do programa, modularização, suas estruturas de dados e saídas de erro. As descrições de cada função implementada estão comentadas diretamente no código fonte.

## Compilação

No terminal, execute: 

```
make
```

## Execução do programa

No terminal, execute: 

```
invmat [-e arquivo_entrada] [-s arquivo_saida] [-r N] -i k
```
Onde:

    -e arquivo_entrada: parâmetro opcional no qual arquivo_entrada é o caminho completo para o arquivo contendo a matriz a ser invertida. Em caso de ausência do parâmetro, a entrada deve ser lida de stdin.
    -s arquivo_saida: parâmetro opcional no qual arquivo_saida é o caminho completo para o arquivo que vai conter a matriz inversa. Em caso de ausência do parâmetro, a saída deve ser impressa em stdout.
    -r N: parâmetro opcional no qual N é a dimensão da matriz de A a ser gerada aleatoreamente de acordo com a função especificada ao final desta página. Com esta opção, a opção -e deve ser ignorada.
    -i k: Número de iterações de refinamento a serem executadas (>0)

Exemplo de entrada via arquivo:

```
./invmat -e teste.in -s teste.out -i 3
```

Exemplo de entrada via parâmetros:

```
./invmat -r 10 -i 3
```

## Modularização e estruturas de dados

- invmat.c
Módulo do programa principal, que contém apenas função main.

São usadas matrizes e vetores alocados dinamicamente para alocar as matrizes e vetores usados no programa, como as matrizes L, U e o vetor LUT, usado para guardar as trocas de linha na fatoração LU e, posteriormente, acessar B.

- sislin.c
Módulo fornecido inicialmente na especificação do trabalho. O módulo sislin.* contém funções para definir um sistema linear. A geração de matrizes aleatórias a serem invertidas, de tamanho NxN, utilizaram as funções alocaSisLin() (tipo pontVet) e iniSisLin(), (tipo generico), especificadas no módulo sislin.* 

A matriz identidade, que é usada na resolução do sistema via Fatoração LU, foi alocada a partir de uma modificação na estrutura do Sistema Linear - e nas respectivas funções que a utilizam -, que transformou SL-b de um vetor em uma matriz de bs (um para cada coluna da matriz identidade):

```
// Estrutura para definiçao de um sistema linear qualquer
typedef struct {
  real_t **A; // coeficientes
  real_t **b; // termos independentes
  unsigned int n; // tamanho do SL
  tipoAloc_t tipoAloc_A; // tipo de alocação usada na matriz de coeficientes
} SistLinear_t;
```
Contém ainda outras novas funções necessárias para implementação e testes:

```
void prnVetor (real_t *v, unsigned int n);
void prnMatriz (real_t **m, unsigned int n);
int copia_matriz(real_t **x, real_t **y, int n);
int copia_vetor(real_t *x, real_t *y, int n);
int copiaSisLin(SistLinear_t *SL, SistLinear_t *SL_copia);
int parseArguments(int argc, char** argv, FILE** fp_in, FILE** fp_out, int *N, int *k, int *flag_e, int *flag_s, int *flag_r, int *flag_i);
real_t** alocaMatriz(int N);
int* alocaeInicilizaVetor(int N);
```
- utils.c
O módulo utils.* contém a definição da função timestamp(), que usada para o cálculo dos tempos de execução. O tempo decorrido é medido pela diferença do "timestamp" medido antes e depois da região de interesse.

- Metodos.c
Esse módulo contém as funções que são o coração do trabalho, ou seja, que fazem de fato os cálculos necessários à inversão das matrizes.

Funções auxiliares:
```
int encontraMaxColunaPivo(double** M, int pivNum, int n);
void criaMatrizIdentidade(real_t **M, int n);
void MultiplicaMQs(double** mA, double** mB, double** mR, int n);
real_t calculaDeterminante(real_t **M, int n);
```

Funções usadas para fatoração LU e cálculo da matriz inversa usando Eliminação de Gauss e refinamento:
```
int FatoracaoLU_PivoParcial(double** L, double** U, int n, int* P, double *tTotal);
void TrocaElementosVetor(int *vet, int i1, int i2);
void trocaLinhasMQ(double **M, int n, int l1, int l2);
void CalculaYFROML(real_t **L, int n, int *LUT, int k, real_t* y);
void CalculaXFROMUY(real_t **U, real_t* y, int n, real_t *x);
int calculaInversa(real_t **L, real_t **U, real_t **I, int *LUT, unsigned int n, double *tTotalY, double *tTotalX);
```

## Problemas

Devido alguns problemas, a equipe não conseguiu terminar o trabalho. Devido a isso, faltou métodos como o refinamento e o calculo de norma de cada iteração (para k iterações, como foi colocado na especificação do trabalho feito pelo professor). Apenas foi calculado a primeira inversa, que por
não ser refinada, tem uma boa diferença da resposta que deveria chegar.

## Pontos positivos

Foi possível explorar bem a parte de erros no que foi implementado, foi feito testes de mesa para acompanha a execução dos algoritmos até o que foi possível, foi possível progredir bastante mesmo devido os problemas, foi modularizado com funções para cada etapa e por pouco que não foi terminado tudo, que foram limpos para uma entrega do que foi possível ser feito.