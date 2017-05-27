import javafx.application.Platform;

public class Tema3 {

    private static double M[][];
    private static double B[];
    private static int m=10;
    private static double e=0.0000000001;
    private static int p=1;
    private static double pi=3.14159265359;

    public static void main(String[] args) {
        //STEP 1: Generate a matrix
        M=generate_matrix();

        //STEP 2: Generate a vector
        B=generate_vector();

        //PRINT MATRIX
        //print_matrix(M);

        //Metoda Jacobi
        //System.out.println("Metoda Jacobi:");
        //metoda_jacobi(M,B);

        //Metoda Gradientului
        System.out.println("Metoda gradientului: ");
        //metoda_gradientului(M,B);

        //Metoda Gauss Seidel
        //System.out.println("Metoda Gauss Seidel:");
        //metoda_gauss_seidel(M,B);

        //Metoda rotatiilor
        System.out.println("Metoda rotatiilor:");
        metoda_rotatiilor(M,B);
    }

    public static void metoda_jacobi(double [][]matrix, double []vector) {
        double maxSum=0;
        for (int i=1;i<=m;i++) {
            double currentSum=0;
            for (int j=1;j<=m;j++) {
                currentSum+=Math.abs(matrix[i][j]);
            }
            if (currentSum>maxSum) {
                maxSum=currentSum;
            }
        }

        double infiniteNorm=maxSum;
        int no=0;
        double []Z=new double[m+1];

        for (int k=1;k<=p;k++) {
            double sigma=(2*k)/(infiniteNorm*(p+1));
            double[][] A=new double[m+1][m+1];
            for (int i=1;i<=m;i++) {
                A[i][i]=1-sigma*matrix[i][i];
            }
            for (int i=1;i<=m;i++) {
                for (int j=1;j<=m;j++) {
                    if (i!=j) {
                        A[i][j]=Math.abs(-sigma*matrix[i][j]);
                    }
                }
            }
            double []bSigma=new double[m+1];
            for (int i=1;i<=m;i++) {
                bSigma[i]=sigma*vector[i];
            }
            double []X=new double[m+1];
            double []Y=new double[m+1];
            int n=0;
            boolean condition=true;

            while (condition) {
                n++;
                for (int i=1;i<=m;i++) {
                   double auxSum=0;
                   for (int j=1;j<=m;j++) {
                       auxSum+=A[i][j]*X[j];
                   }
                   Y[i]=auxSum+bSigma[i];
                }
                double auxSum=0;
                for (int i=1;i<=m;i++) {
                    for (int j=1;j<=m;j++) {
                        auxSum+=matrix[i][j]*(Y[i]-X[i])*(Y[j]-X[j]);
                    }
                }
                for (int i=1;i<=m;i++) {
                    X[i]=Y[i];
                }
                condition=(Math.sqrt(auxSum)>=e);
            }
            if (k==1 || n<no) {
                no=n;
                Z=X;
            }
        }
        System.out.println("Numarul de iteratii="+no);
        for (int i=1;i<=m;i++) {
            System.out.printf("%.4f ",Z[i]);
        }
        System.out.println();
        System.out.println();
    }

    public static void metoda_gauss_seidel(double [][]matrix, double []vector) {
        int no=0;
        double []Z=new double[m+1];
        for (int k=1;k<=p;k++) {
            double sigma=(2*k)/(p+1);
            double X[]=new double[m+1];
            double Y[]=new double[m+1];
            int n=0;
            boolean condition=true;
            while (condition) {
                n++;
                for (int i=1;i<=m;i++) {
                    double auxSum1=0,auxSum2=0;
                    for (int j=1;j<=i-1;j++) {
                        auxSum1+=matrix[i][j]*Y[j];
                    }
                    for (int j=i+1;j<=m;j++) {
                        auxSum2+=matrix[i][j]*X[j];
                    }
                    Y[i]=(1-sigma)*X[i]+(sigma/matrix[i][i])*(vector[i]-auxSum1-auxSum2);
                }
                double auxSum=0;
                for (int i=1;i<=m;i++) {
                    for (int j=1;j<=m;j++) {
                        auxSum+=matrix[i][j]*(Y[i]-X[i])*(Y[j]-X[j]);
                    }
                }
                for (int i=1;i<=m;i++) {
                    X[i] = Y[i];
                }
                condition=(Math.sqrt(auxSum)>e);
            }
            if (k==1 || n<no) {
                no=n;
                Z=X;
            }
        }
        System.out.println("Numarul de iteratii="+no);
        for (int i=1;i<=m;i++) {
            System.out.printf("%.4f ",Z[i]);
        }
        System.out.println();
        System.out.println();
    }

    public static void metoda_gradientului(double [][]matrix, double []vector) {
        double []X0=new double[m+1];
        double []R0=multiplicationMatrix(matrix,X0);
        for (int i=1;i<=m;i++) {
            R0[i]=vector[i]-R0[i];
        }
        double []V0=R0;
        int n=0;
        boolean condition=true;
        while (condition) {
            double auxSum1=0;
            for (int i=1;i<=m;i++) {
                auxSum1+=R0[i]*R0[i];
            }
            double []Av=multiplicationMatrix(matrix,V0);
            double auxSum2=0;
            for (int i=1;i<=m;i++) {
                auxSum2+=Av[i]*V0[i];
            }
            double a=auxSum1/auxSum2;
            double []X1=new double[m+1];
            double []R1=new double[m+1];
            for (int i=1;i<=m;i++) {
                X1[i]=X0[i]+a*V0[i];
                R1[i]=R0[i]-a*Av[i];
            }
            auxSum1=0;
            for (int i=1;i<=m;i++) {
                auxSum1+=R1[i]*R1[i];
            }
            auxSum2=0;
            for (int i=1;i<=m;i++) {
                auxSum2+=R0[i]*R0[i];
            }
            double c=auxSum1/auxSum2;
            double []V1=new double[m+1];
            for (int i=1;i<=m;i++) {
                V1[i]=R1[i]+c*V0[i];
            }
            n++;
            double auxSum=0;
            double X1mX0[]=new double[m+1];
            for (int i=1;i<=m;i++) {
                X1mX0[i]=X1[i]-X0[i];
            }
            for (int i=1;i<=m;i++) {
                auxSum+=X1mX0[i]*X1mX0[i];
            }
            for (int i=1;i<=m;i++) {
                X0[i] = X1[i];
                V0[i] = V1[i];
                R0[i] = R1[i];
            }
            condition=(Math.sqrt(auxSum)>=e);
        }
        System.out.println("Numarul de iteratii="+n);
        for (int i=1;i<=m;i++) {
            System.out.printf("%.4f ",X0[i]);
        }
        System.out.println();
        System.out.println();
    }

    public static void metoda_rotatiilor(double [][]matrix, double []vector) {
        double [][]X=matrix;
        double [][]Y=new double[m+1][m+1];
        int n=0;
        boolean condition=true;
        while (condition) {
            int p=1,q=1;
            double maxim=Math.abs(X[1][1]);
            for (int i=1;i<=m;i++) {
                for (int j=1;j<=m;j++) {
                    if (i<j) {
                        if (Math.abs(X[i][j])>maxim) {
                            p=i;
                            q=j;
                            maxim=Math.abs(X[i][j]);
                        }
                    }
                }
            }
            double tetha;
            if (X[p][p]==X[q][q]) {
                tetha=pi/4;
            } else {
                tetha=1/2*Math.atan((2*X[p][q])/X[p][p])-X[q][q];
            }
            double c=Math.cos(tetha);
            double s=Math.sin(tetha);
            for (int i=1;i<=m;i++) {
                for (int j=1;j<=m;j++) {
                    if (i!=p && i!=q && j!=p && j!=q) {
                        Y[i][j]=X[i][j];
                    }
                }
            }
            for (int j=1;j<=m;j++) {
                if (j!=p && j!=q) {
                    Y[p][j]=c*X[p][j]+s*X[q][j];
                    Y[j][p]=c*X[p][j]+s*X[q][j];
                    Y[q][j]=c*X[q][j]-s*X[p][j];
                    Y[j][q]=c*X[q][j]-s*X[p][j];
                }
            }
            Y[p][q]=0;
            Y[q][p]=0;
            Y[p][p]=c*c*X[p][p]+2*c*s*X[p][q]+s*s*X[q][q];
            Y[q][q]=s*s*X[p][p]-2*c*s*X[p][q]+c*c*X[q][q];
            double auxSum=0;
            for (int i=1;i<=m;i++) {
                for (int j=1;j<=m;j++) {
                    if (i!=j) {
                        auxSum+=Y[i][j]*Y[i][j];
                    }
                }
            }
            double modul=Math.sqrt(auxSum);
            for (int i=1;i<=m;i++) {
                for (int j=1;j<=m;j++) {
                    X[i][j]=Y[i][j];
                }
            }

            n++;
            condition=(modul>=e) && n<500;
        }
        System.out.println("Valorile proprii ale matricei: "+n);
        for (int i=1;i<=m;i++) {
            for (int j=1;j<=m;j++) {
                System.out.print(X[i][j]+" ");
            }
            System.out.println();
        }
        System.out.println();
        System.out.println();
    }

    private static double[][] generate_matrix() {
        double matrix[][]=new double[m+1][m+1];
        for (int i=1;i<=m;i++) {
            matrix[i][i]=2+1.0/(m*m);
        }
        for (int i=1;i<=m-1;i++) {
            matrix[i][i+1]=matrix[i+1][i]=-1;
        }
        return matrix;
    }

    private static double[] generate_vector() {
        double temp[]=new double[m+1];
        for (int i=1;i<=m;i++) {
            temp[i]=1.0/(m*m);
        }
        return temp;
    }

    private static void print_matrix(double[][] matrix) {
        for (int i=1;i<matrix.length;i++) {
            for (int j=1;j<matrix.length;j++) {
                System.out.printf("%.2f ",matrix[i][j]);
            }
            System.out.println();
        }
        System.out.println();
    }

    private static double[] multiplicationMatrix(double[][] m1,double[] m2) {
        double[] temp=new double[m1.length];
        for (int i = 1; i < m1.length; i++) {
            for (int j = 1; j < m1.length; j++) {
                temp[i] += m1[i][j] * m2[i];
            }
        }
        return temp;
    }
}
