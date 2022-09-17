# Trabalho 1 - Introdução à Computação Científica: inversão de matriz usando fatoração LU

## Equipe

- João Victor Kenji Tamaki - GRR20182965 
- Marisa Sel Franco - GRR20186556

## Descrição

-> Melhoria de Desempenho

Você deve alterar o código do primeiro trabalho (v1) de forma a obter uma melhora no desempenho (v2) em duas partes:

Na operação de resolução do sistema linear triangular (op1): LUx=b -> Ly=b -> Ux=y;
Na operação de cálculo do resíduo (op2): R = I - A * inv(A);
As alterações devem ser explicadas no relatório a ser entregue, justificando as razões pelas quais você efetuou cada alteração.

MUITO IMPORTANTE: você deve ser explicar por que motivo suas modificações levaram a um aumento de desempenho.


-> Análise de Desempenho

Uma vez alterado o código, você deve comparar o desempenho das duas versões. Estas análises devem ser descritas, sob a forma de um trabalho acadêmico, e entregues em formato PDF.

É imprescindível que sejam respeitadas as seguintes condições:

Ambos códigos devem ser compilados com GCC e as opções: -O3 -mavx -march=native;
Os códigos devem ser compilados na mesma máquina utilizada para os testes.
Os testes devem utilizar os mesmos parâmetros e em igualdade de condições;
Ambos códigos devem ser instrumentados com a biblioteca do LIKWID para que se possa separar o cálculo da op1 do cálculo da op2, ou seja, você deve criar um marcador LIKWID para a parte da op1 e outro para a parte da op2;
Você pode escolher um computador de sua preferência, desde que possua os contadores Likwid especificados. Não utilize as servidoras de processamento do DInf que tenham uso compartilhado. Elas podem ser máquinas virtuais e o compartilhamento impede medidas de desempenho. Em caso de dúvida, consulte o professor.
Você deve apresentar a arquitetura do processador utilizado nos testes no seu texto. Estas informações podem ser obtidas com o comando

         likwid-topology -g -c
 Para comparar o desempenho dos códigos, você deve efetuar uma série de testes.

Cada teste deve ser reportado sob a forma de um gráfico de linhas, com linhas em cores distintas para os resultados de cada versão (v1 e v2);
No eixo das abcissas (eixo x) os gráficos representam o tamanho da matriz (N).
Cada teste deve ser executado para os seguintes tamanhos de matriz:  N={32, 33, 64, 65, 128, 129, 256, 257, 512, 1000, 2000, 4000 6000 10000};
O número de iterações em todos os testes, para todos tamanhos de matrizes, deve ser -i 10;
As matrizes devem ser geradas aleatoreamente com a função definida no primeiro trabalho;
Cada teste deve apresentar em linhas distintas os valores para o cálculo de cada operações (op1 e op2). Assim, os gráficos terão sempre 4 linhas, duas para a v1 e duas para a v2;
Cada gráfico deve ser explicado e você deve demonstrar que consegue entender o que está reportado nele;
Os gráficos devem ser apresentados com o eixo das ordenadas (eixo y) em escala logarítmica.
 

Os seguintes testes devem ser executados (um gráfico para cada teste):

Teste de tempo: mostra o tempo médio do cálculo da op1 e o tempo médio do cálculo da op2 (utilize a função "timestamp()" para medir o tempo);
Banda de Memória: utilizar o grupo MEM do LIKWID, e apresentar o resultado de "Memory bandwidth [MBytes/s]"; Caso não tenha o grupo MEM, utilize o grupo L3.
Cache miss L1: utilizar o grupo CACHE ou L1CACHE do LIKWID, e apresentar o resultado de "data cache miss ratio".
Caso não tenha o cache miss da L1, utilize o cache miss da L2 (grupo L2CACHE)
Operações aritméticas: utilizar o grupo FLOPS_DP ou FLOPS_AVX do LIKWID, e apresentar o resultado de "MFLOP/s"


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

## Otimizações realizadas

Já efetuadas no T1:

-> Escolha do Layout Structure-of-arrays: 
-> Acesso contínuo de matrizes (linha a linha) para evitar stride
-> Lookup Table para o pivoteamento parcial, auxiliando no desempenho do cálculo da inversa.

Efetuadas para o T2:

-> Uso eficiente do pipeline: Sem dependência dos dados e argumentos passados como Restrict para indicar que aliasing não foi usado;
-> Loop unroll & Jam para otimizar o pipeline, maximizar o uso da cache e a quantidade de operações em cada loop, diminuindo cache miss (menos acesso a memória) e o tempo de execução, pois aproveita melhor o pipeline;
-> Matrizes alocadas contiguamente na memória (com alligned_alloc) e com padding e cuidados de tamanho das dimensões para evitar cache trashing;
-> Economizando operações aritméticas por meio do cálculo das linhas que serão usadas em cada iteração de linha, evitando realizar esse cálculo durante toda a coluna e reduzindo o tempo de cada iteração do loop;
-> Foram criados scripts para automatizar os testes;
-> Utilizada as flags para compilação necessárias para a otimização e vetorização e utilização do LIKWID;
-> Não foi utilizado inline para não haver mais chances de sobrecarregar os registradores, devido já estar sendo utilizado o Loop Unroll & Jam;
-> Organização, acesso e utilização das estruturas de dados para melhorar o uso SIMD, seguindo as instruções dadas em aula;
-> Código auxiliar em Python para automatizar a geração de Gráficos.

## Problemas

Devido ao acesso remoto, era possível cair a rede por algum tempo ou algum outro usuário entrar e competir com o uso dos recursos do computador. Devido a equipe só ter um integrante ativo, não foi possível realizar mais testes e tentar otimizar mais o código.

## Conclusão

Neste trabalho foi possível testar e praticar bastante o uso das técnicas de otimização e foi atingido um resultado muito promissor. Mesmo que o tempo para execução tenha sido apertado, foi possível terminar a implementação do primeiro trabalho e fazer muitas alterações para otimização do código.