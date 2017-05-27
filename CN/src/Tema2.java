
public class Tema2 {

    private static int M[][];
    private static int[] B;
    private static int n=5;
    private static int p=10;

    public static void main(String[] args) {
        //STEP 1 --- GENERATE A MATRIX
        M=generate_matrix(n,p);

        System.out.println("Matricea generata:");
        print_matrix(M);

        //STEP2 --- GENERATE VECTOR B
        B=generate_vector(n,1);

        //Metoda LU
        // @INFO: https://en.wikipedia.org/wiki/LU_decomposition
        System.out.println("Metoda LU:");
        metoda_LU(n,M,B);

        //Metoda Cholesky
        // @INFO: https://rosettacode.org/wiki/Cholesky_decomposition#Java
        System.out.println("Metoda Cholesky:");
        int MT[][]=transposeMatrix(M);
        int M2[][]=multiplicationMatrix(MT,M);
        int B2[]=multiplicationMatrix2(MT,B);
        metoda_Cholesky(n,M2,B2);

        //Metoda QR
        // @INFO: https://ro.wikipedia.org/wiki/Descompunerea_QR
        System.out.println("Metoda QR:");
        metoda_QR(n,M,B);

    }

    public static void metoda_LU(int n, int [][]matrix, int [] vector) {
        double [][]L=new double[n+1][n+1];
        double [][]U=new double[n+1][n+1];

        for (int i=1;i<=n;i++) {

            //CREATE L
            for (int j = i; j <= n; j++) {
                double auxSum = 0;
                for (int k = 1; k <= i - 1; k++) {
                    auxSum += L[j][k] * U[k][i];
                }
                L[j][i] = matrix[j][i] - auxSum;
            }

            //CREATE U
            for (int j=i+1;j<=n;j++) {
                double auxSum = 0;
                for (int k=1;k<=i-1;k++) {
                    auxSum+=L[i][k]*U[k][j];
                }
                U[i][j]=(matrix[i][j]-auxSum)/L[i][i];
            }
        }

        double []x=new double[n+1];
        for (int i=1;i<=n;i++) {
            int auxSum=0;
            for (int j=1;j<=i-1;j++) {
                auxSum+=L[i][j]*x[j];
            }
            x[i]=(vector[i]-auxSum)/L[i][i];
        }

        double []y=new double[n+1];
        for (int i=n;i>=1;i--) {
            int auxSum=0;
            for (int j=i+1;j<=n;j++) {
                auxSum+=U[i][j]*y[j];
            }
            y[i]=x[i]-auxSum;
        }

        for (int i=1;i<=n;i++) {
            System.out.printf("%.0f ",y[i]);
        }
        System.out.println();
        System.out.println();
    }

    public static void metoda_Cholesky(int n, int[][] matrix, int[] vector) {
        double [][]L=new double[n+1][n+1];
        for (int i=1;i<=n;i++) {
            double auxSum=0;
            for (int j=1;j<=i-1;j++) {
                auxSum+=L[i][j]*L[i][j];
            }

            if (matrix[i][i]-auxSum<=0) {
                System.out.println("Conditia nu este indeplinita!");
                System.exit(-1);
            } else {
                L[i][i]=Math.sqrt(matrix[i][i]-auxSum);
            }
            for (int j=i+1;j<=n;j++) {
                auxSum=0;
                for (int k=1;k<=i-1;k++) {
                    auxSum+=L[j][k]*L[i][k];
                }
                L[j][i]=(matrix[j][i]-auxSum)/L[i][i];
            }
        }

        double []x=new double[n+1];
        for (int i=1;i<=n;i++) {
            double auxSum=0;
            for (int j=1;j<=i-1;j++) {
                auxSum+=L[i][j]*x[j];
            }
            x[i]=(vector[i]-auxSum)/L[i][i];
        }

        double []y=new double[n+1];
        for (int i=n;i>=1;i--) {
            double auxSum=0;
            for (int j=i+1;j<=n;j++) {
                auxSum+=L[j][i]*y[j];
            }
            y[i]=(x[i]-auxSum)/L[i][i];
        }

        for (int i=1;i<=n;i++) {
            System.out.printf("%.0f ",y[i]);
        }
        System.out.println();
        System.out.println();
    }

    public static void metoda_QR(int n, int[][] matrix, int[] vector) {
        double [][]Q=new double[n+1][n+1];
        double [][]R=new double[n+1][n+1];
        double auxSum=0;
        for (int i=1;i<=n;i++) {
            auxSum+=matrix[i][1]*matrix[i][1];
        }
        R[1][1]=Math.sqrt(auxSum);
        if (R[1][1]==0) {
            System.exit(-1);
        }

        for (int i=1;i<=n;i++) {
            Q[i][1]=matrix[i][1]/R[1][1];
        }

        for (int k=1;k<=n;k++) {
            for (int j=1;j<=k-1;j++) {
                auxSum=0;
                for (int i=1;i<=n;i++) {
                    auxSum+=matrix[i][k]*Q[i][j];
                }
                R[j][k]=auxSum;
            }

            double auxSum1=0,auxSum2=0;
            for (int i=1;i<=n;i++) {
                auxSum1+=matrix[i][k]*matrix[i][k];
            }
            for (int i=1;i<=k-1;i++) {
                auxSum2+=R[i][k]*R[i][k];
            }
            R[k][k]=Math.sqrt(auxSum1-auxSum2);
            for (int i=1;i<=n;i++) {
                auxSum=0;
                for (int j=1;j<=k-1;j++) {
                    auxSum+=R[j][k]*Q[i][j];
                }
                Q[i][k]=(matrix[i][k]-auxSum)/R[k][k];
            }
        }

        double []x=new double[n+1];
        for (int i=1;i<=n;i++) {
             auxSum=0;
            for (int j=1;j<=n;j++) {
                auxSum+=Q[j][i]*vector[j];
            }
            x[i]=auxSum;
        }

        double []y=new double[n+1];
        for (int i=n;i>=1;i--) {
             auxSum=0;
            for (int j=i+1;j<=n;j++) {
                auxSum+=R[i][j]*y[j];
            }
            y[i]=(x[i]-auxSum)/R[i][i];
        }

        for (int i=1;i<=n;i++) {
            System.out.printf("%.0f ",y[i]);
        }
        System.out.println();
        System.out.println();
    }

    private static int[] generate_vector(int n, int value) {
        int [] vector=new int[n+1];
        for (int i=1;i<=n;i++) {
            vector[i]=value;
        }
        return vector;
    }

    private static int[][] generate_matrix(int n, int p) {
        int matrix[][]=new int[n+1][n+1];
        for (int i=1;i<=n;i++) {
            for (int j=1;j<=n;j++) {
                matrix[i][j]=binomial(p+j-1,i-1);
            }
        }
        return matrix;
    }

    private static void print_matrix(int matrix[][]) {
        for (int i=1;i<matrix.length;i++) {
            for (int j=1;j<matrix.length;j++) {
                System.out.print(matrix[i][j]+" ");
            }
            System.out.println();
        }
        System.out.println();
    }

    private static int binomial(int n, int k) {
        if (k>n-k)
            k=n-k;

        int b=1;
        for (int i=1, m=n; i<=k; i++, m--)
            b=b*m/i;
        return b;
    }

    private static int[][] transposeMatrix(int [][] m){
        int[][] temp = new int[m.length][m.length];
        for (int i = 1; i < m.length; i++)
            for (int j = 1; j < m.length; j++)
                temp[j][i] = m[i][j];
        return temp;
    }

    private static int[][] multiplicationMatrix(int[][] m1,int[][] m2) {
        int[][] temp = new int[m1.length][m1.length];
        for (int i = 1; i < m1.length; i++) {
            for (int j = 1; j < m1.length; j++) {
                for (int k = 1; k < m1.length; k++) {
                    temp[i][j] += m1[i][k] * m2[k][j];
                }
            }
        }
        return temp;
    }

    private static int[] multiplicationMatrix2(int[][] m1,int[] m2) {
        int[] temp=new int[m1.length];
        for (int i = 1; i < m1.length; i++) {
            for (int j = 1; j < m1.length; j++) {
                    temp[i] += m1[i][j] * m2[i];
            }
        }
        return temp;
    }
}
