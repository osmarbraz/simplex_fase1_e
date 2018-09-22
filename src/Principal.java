
/*************************************************************************
 * Considera somente problemas que possam ser resolvidos com a fase 1 
 * do Algoritmo Simplex.
 * 
 * Dado uma matriz A m por n , um vetor B de tamanho M e um
 * vetor c de tamanho n, resolve um problema de programação linear
 * { max cx : Ax <= b, x >= 0 }. Assume que b >= 0 e que x = 0 é a solução
 * básica factivel.
 *
 * Cria um tableau simplex (m+1)por (n+m+1) com a coluna RHS n+m, a função
 * objetivo na linha m, e as variáveis de folga nas colunas m a n+m+1
 * 
 * M = {1,2,..., m}, o conjunto dos índices das restrições do problema
 * N = {1,2,..., n}, o conjunto dos índices das variáveis do problema
 * z é a função objetivo;
 * xj são as variáveis de decisão, principais ou controláveis, j = 1, 2,..., n;
 * aij é a constante ou coeficiente da i-ésima restrição da j-ésima variável, i = 1, 2,...,m, j = 1, 2,...,n;
 * bi é o termo independente ou quantidade de recursos disponíveis da i-ésima restrição, i = 1, 2,...,m;
 * cj é a constante ou coeficiente da j-ésima variável da função objetivo, j = 1, 2,...,n
 * 
 *************************************************************************/

public class Principal {
    
    private static double[][] a;   // tableau
    private static int m;          // número de restrições
    private static int n;          // número de variáveis não básicas
    private static int iteracoes ; // número de iterações
    private static int[] base;     // base[i] = variável básica correspondente a linha i
                                    // Necessário somente para imprimir a solução.
    
    /**
     * Monta tableau original.
     * 
     * @param A Matriz das restrições.
     * @param b Vetor das constantes de A.
     * @param c Coeficientes da função objetivo.
     */    
    public static void montaTableau(double[][] A, double[] b, double[] c){
        //Quantidade de restrições do problema
        m = b.length;
        //Quantidade de variáveis do problema
        n = c.length;
        //Uma linha e coluna adicional para função objetivo
        a = new double[m + 1][n + m + 1];
        //Copia as restrições para a matriz A para a
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                a[i][j] = A[i][j];
            }
        }
        //Adiciona 1 para a diagonal principal das variaveis não básicas
        //Considera a partir de n+1 onde n é a quantidade de variáveis 
        //básicas até m que é a quantidade de variáveis não básicas
        for (int i = 0; i < m; i++) {
            a[i][n + i] = 1.0;
        }
        //Adiciona a constante da função objetivo para a linha m
        for (int j = 0; j < n; j++) {
            //Inverte o sinal do coeficiente c
            a[m][j] = -c[j];
        }
        
        //Adiciona constante das restrições a coluna n+m
        for (int i = 0; i < m; i++) {
            a[i][n + m] = b[i];
        }

        //Vetor auxiliar para imprimir a solução das variáveis básicas
        base = new int[m];
        for (int i = 0; i < m; i++) {
            base[i] = n + i;            
        }
    }
    
    /** 
     * Roda o algoritmo simplex apartir de uma solução inicial BFS.
     * 
     * @param A Matriz das restrições.
     * @param b Vetor das constantes de A.
     * @param c Coeficientes da função objetivo.
     */
    public static void solver(double[][] A, double[] b, double[] c) {
        //Passo 1 - Monta o Tableau
        montaTableau(A,b,c);
        
        //Conta as iterações
        iteracoes = 0;
        // Passo 2 - Verifica condições de otimalidade
        while (testarOtimalidade(c)==true) {
            
            //Passo 2.1 Regra de Entrada 
            //encontra a coluna que deve entrar na base
            int q = encontraVariavelEntrada();
            
            //Passo 2.2 Regra de Saída
            // encontra a linha que deve sair considerando a variável de entrada
            int p = encontraVariavelSaida(q);
            
            //Passo 3 Recalcula a base
            recalculaBase(p, q);
            
            // atualiza a base guardando a variável que entrou na base.
            if (p!=-1){
                base[p] = q;
            }
            
            //Incrementa o contador de iterações
            iteracoes = iteracoes + 1;            
        }                     
    }

    /**
     * Indice de uma coluna não básica com o maior custo.
     * 
     * @return Coluna da variável de entrada ou -1.
     */    
    private static int encontraVariavelEntrada() {
       double maior = 0;
        int colunaMaior = -1;
        //Procura somente nas variáveis não base
        for (int j = 0; j < n; j++) {            
            //Usa o primeiro elemento como maior para inicializar p
            if (colunaMaior==-1){
                maior =  Math.abs(a[m][j]);
                colunaMaior = j;
            } else {
                //Maior maior absoluto
                if (Math.abs(a[m][j]) > maior) {
                    maior =  Math.abs(a[m][j]);
                    colunaMaior = j;
                }
            }
        }        
        return colunaMaior;  
    }
  
    /**
     * Encontra linha p usando a regra da menor razão.(-1 se não houver a linha).
     * 
     * @param k Coluna do pivô.
     * @return A linha da variável de saída ou -1.
     */      
    private static int encontraVariavelSaida(int k) {
        double menor = 0;
        int linhaMenor = -1;        
        // Posição do pivô deve ser diferente de -1
        if (k != -1){
            for (int i = 0; i < m; i++) {
                //aik deve ser positivo
                if (a[i][k] > 0) {
                    //Usa o primeiro elemento maior que 0 como menor para inicializar p
                    if (linhaMenor==-1){
                        linhaMenor = i;
                        menor = a[i][m + n] / a[i][k];
                    } else {
                        if ((a[i][m + n] / a[i][k]) < menor) {
                            linhaMenor = i;
                            menor =  a[i][m + n] / a[i][k];
                        }
                    }               
                } 
            }
        }
        return linhaMenor;
    }

    /**
     * Calcula o pivô com a entrada(p, q) usando o método de eliminação de Gauss-Jordan.
     * 
     * @param p Linha do pivô.
     * @param q Coluna do pivô.
     */
    private static void recalculaBase(int p, int q) {
        //Calcula somente se existir um pivô        
        // p e q devem ser diferente de -1
        if ((p!=-1) && (q != -1)){    
            
            //Calcula a linha do pivô
            //Guarda o valor do pivô
            double valorPivo = a[p][q];
            //Percorre a coluna
            for (int j = 0; j <= n + m; j++) {
                //Recalcula o elemento da linha a[p][j]
                a[p][j] = a[p][j]/valorPivo;
            }
            
            // Recalcula as outras linhas da matriz a menos a linha p            
            //Percorre as linhas
            for (int i = 0; i <= m; i++) {
                //Menos para linha do pivô
                if (i != p) {      
                    //Guarda o coeficiente da coluna do pivô
                    double coeficienteColunaPivo = a[i][q];
                    //Percorre a coluna
                    for (int j = 0; j <= n + m; j++) {
                        //Recalcula o elemento
                        a[i][j] = a[i][j] - (coeficienteColunaPivo * a[p][j]);
                    }
                }
            }
        }
    }

    /**
     * Retorna o valor ótimo da função objetivo.
     * 
     * @return Valor ótimo da função objetivo.
     */
    public static double getValor() {
        //Última linha e coluna da matriz a
        return a[m][n + m];
    }

     /**
     * Testa a otimalidade verificando se é a solução básica fáctivel.
     * 
     * Procura algum valor negativo em c(linha m).
     * 
     * Como os coeficientes de x1 e x2  são negativos na linha m, 
     *  a SBF(Solução Básica Factível) atual não é ótima, pois um 
     *  incremento positivo em x1  ou x2  resultará em SBF adjacente 
     *  melhor do que a SBF atual.
    */
    private static boolean testarOtimalidade(double[] c) {
        int k = 0;
        boolean temNegativo = false;
        //Se existir um elemento negativo interrompe o laço
        while ((k < c.length) && (temNegativo==false)){
            // verifica se a[m][k] < 0
            // m é a última linha da matriz a
            if (a[m][k] < 0.0){
                temNegativo = true;                
            }
            k = k + 1;
        }        
        return temNegativo;
    }

    /**
     * Mostra o tableau.
     */
    public static void mostraTableau() {
        System.out.println("m = " + m + " e n = " + n);        
        //Cabeçalho do tableau
        System.out.printf("\t\t");
        for (int i = 0; i < n + m; i++) {
            System.out.printf("   x[%d]\t\t",i);
        }
        System.out.printf("      b");
        System.out.println();

        //Matriz a        
        for (int i = 0; i <= m; i++) {
            if (i>=m){
                System.out.printf("z =\t");
            } else {
                System.out.printf("x["+(n+i)+"] =\t");
            }
            for (int j = 0; j <= n + m; j++) {
                System.out.printf("\t%7.2f ", a[i][j]);
            }
            System.out.println();
        }    
    }    
    
    /**
     *  retorna o vetor da solução primal.
     * 
     * @return um vetor com a solução
     */
    public static double[] getPrimal() {
        double[] x = new double[n];
        for (int i = 0; i < m; i++) {
            if (base[i] < n) {
                x[base[i]] = a[i][n + m];
            }
        }
        return x;
    }
    
     /**
     * Imprime a solução do problema.
     */
    public static void imprimeSolucao(){
        //Função objetivo        
        System.out.println("A solução foi encontrada com " + iteracoes + " iterações");
        System.out.println("   z = " + getValor());
        //Variáveis básicas
        double[] x = getPrimal();
        for (int i = 0; i < x.length; i++) {
            System.out.println("x[" + (i+1) + "] = " + x[i]);
        }                
    }
    
    /**
     * Executa o algoritmo solver.
     * 
     * @param A Matriz das restrições.
     * @param b Vetor das constantes de A.
     * @param c Coeficientes da função objetivo.
     */
    public static void teste(double[][] A, double[] b, double[] c) {        
        //Monta o tableau
        montaTableau(A, b, c);
        
        System.out.println("\n>>>> Tableau Inicial");                
        mostraTableau();
        
        //Executa o solver
        solver(A, b, c);
        
        System.out.println("\n>>>> Tableau Final");                
        mostraTableau();
        
        //Imprime a solução
        imprimeSolucao();
    }

    /**
     * Caso 1
     * 
     * MAX 3X1 + 5X2
     * 
     * SUBJECT TO
     * 1X1 + 0X2 <=4
     * 0X1 + 2X2 <=12
     * 3X1 + 2X2 <=18
     * X1 >= 0
     * X2 >= 0
     * END
     * 
     *  3   5   0   0   0   0
     *  1   0   1   0   0   4
     *  0   2   0   1   0   12
     *  3   2   0   0   1   18
     */
    public static void testeCaso1() {
        //Matriz A das equações de restrições
        double[][] A = {
            {1, 0},
            {0, 2},
            {3, 2}
        };        
        //Coeficientes da função objetivo
        double[] c = {3, 5};
        //Constante das equações de A 
        double[] b = {4, 12, 18};
        //Executa o teste para o problema
        teste(A, b, c);
    }
    
    
    /**
     * Caso 2 
     * 
     * Max 1000x1 + 1800x2
     *
     * SUBJECT TO 
     * 20x1 + 30x2 <=1200 
     * x1 <=40 
     * x2 <=30
     * x1 >= 0 x2 >= 0 
     * END
     * 
     *  1000   1800   0    0    0      0
     *    20     30   1    0    0   1200
     *     1      0   0    1    0     40
     *     0      1   0    0    1     30
     */
    public static void testeCaso2() {
        //Matriz A das equações de restrições
        double[][] A = {
            {20, 30},
            {1, 0},
            {0, 1,}
        };        
        //Constante das equações de A 
        double[] b = {1200, 40, 30};
        //Coeficientes da função objetivo
        double[] c = {1000, 1800};
        //Executa o teste para o problema
        teste(A, b, c);
    }
   
    /**
     * Programa principal.
     * 
     * @param args 
     */
    public static void main(String[] args) {
        
        //Executa o teste do caso 1
        testeCaso1();        
        System.out.println("--------------------------------");
        
        //Executa o teste do caso 2
        testeCaso2();        
        System.out.println("--------------------------------");        
    }
}
