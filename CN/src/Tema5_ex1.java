/*
 * @author Cristina Ciort
 */
import java.util.ArrayList;
import java.util.Scanner;

class InitMatrice 
{

    double[][] matrJacobi ;
    double[] b1;
    double[] x1;
    int n;
    double precizie;

    public InitMatrice(int n, double precizie, double[] a, double[] b, double[] c,double w, double z)
    {
    	this.n = n;
    	this.precizie = precizie;
    	this.matrJacobi = new double[n+1][n+1];
    	this.b1 = new double[n+1];
    	this.x1 = new double[n+1];
    	
    	for(int i = 1; i <= n; i++)
    	{
            for(int j = 1; j <= n; j++)
            {
                matrJacobi[i][j] = 0.0;
            }
    	}
    
    	for(int i=2;i<n;i++)
    	{
            matrJacobi[i][i] = 2.0;
    	}
    	
    	matrJacobi[1][1] = 1.0;
    	matrJacobi[n][n] = 1.0;
   
    	
    	for (int i = 2; i <= n-1 ; i++) 
    	{
            matrJacobi[i][i+1] = b[i];
            matrJacobi[i+1][i] = a[i];
    	}
    	
        b1[1] = w;
        b1[n] = z;
        
        for(int i=2;i<n;i++)
            b1[i] = c[i-1];
    }

    public static double normaInfinit(double[][] A ,int n)
    {
    	double niMAX = Double.MIN_VALUE;
    	for(int i = 1; i <= n; i++)
    	{
            double ni = 0.0;
            for(int j = 1; j <= n; j++)
            {
                ni += Math.abs(A[i][j]);
            }
            niMAX = Math.max(niMAX, ni);
    	}
    	return niMAX;
    }

    public double[] Jacobi(int p)
    {
    	double SIGMA;
    	double[][] Bsigma = new double[n+1][n+1];
    	double[] bsigma = new double[n+1];
    	int moptim = 0;
    	double ni = normaInfinit(matrJacobi,n);
    	double l;
    
    	l = 2.0 / ni;
    	for (int k = 1; k <= p-1 ; k++) 
    	{  
            SIGMA = ((double)l / (double)p) * ((double)k);

            for (int i = 1; i <= n; i++)
            {
                for (int j = 1; j <= n; j++)
                {
                    if (i == j)
                    {
                        Bsigma[i][j] = 1.0 - SIGMA * matrJacobi[i][i];
                    } 
                    else
                    {
                        Bsigma[i][j] = -1.0 * SIGMA * matrJacobi[i][j];
                    }
                }
            }

            for (int i = 1; i <= n; i++)
            {
                bsigma[i] = SIGMA * b1[i];
            }

            double norma;
            int m = 0;
            double[] x0 = new double[n+1];
            double[] y  = new double[n+1];
            double s;

            for (int i = 1; i <= n; i++)
            {
                x0[i] = 0.0;
            }

            do 
            {
                m++;
                for (int i = 1; i <= n; i++)
                {
                    s = 0;
                    for (int j = 1; j <= n; j++)
                    {
                        s += (Bsigma[i][j] * x0[j]);
                    }
                y[i] = s + bsigma[i];
                }

                norma = 0;
                for (int i = 1; i <= n; i++)
                {
                    for (int j = 1; j <= n; j++)
                    {
                        norma += matrJacobi[i][j] * (y[j] - x0[j]) * (y[i] - x0[i]);
                    }
                }

                    System.arraycopy(y, 0, x0, 0, n+1);

            } while (norma > precizie);


            if (k == 1)
            {
                moptim = m;
                System.arraycopy(x0, 0, x1, 0, n+1);
            }

            if (k > 1)
            {
                if (m < moptim)
                {
                    moptim = m;
                    System.arraycopy(x0, 0, x1, 0, n+1);
                }
            }
    }
    	
	  return x1;

    }
}


class Puncte implements Comparable<Puncte> {

    double X, Y;

    public Puncte(double left, double right) {
        this.X = left;
        this.Y = right;
    }

    @Override
    public String toString() {
        return ("(" + X + ", " + Y + ")");
    }

    @Override
    public int compareTo(Puncte P) {
        if (Math.abs(this.X - P.X) > 0.001 || Math.abs(this.Y - P.Y) > 0.001 )
            return 1;
        return 0;
 
    }
}

class Interpolare_spline 
{

    private int n, m;
    private double[] V, X, Y, H, A, B, C, Alpha, Beta;
    private double w, z, a, b;
    private ArrayList<Puncte> Solutie;
    
    public static double f(double x) {
        return x * Math.sin(x) + (x * x + 4) * Math.exp(x) - Math.cos(x);
    }
    
    public static double f_deriv(double x){
        return Math.sin(x) + x * Math.cos(x) + 2 * x * Math.exp(x) + (x * x + 4) * Math.exp(x) + Math.sin(x);
    }

    public Interpolare_spline(double[] pct, double a1, double b1, int m1, int n)
    {
        this.n     = n;
        this.X     = new double[n+1];
        this.V     = new double[n+1];
        this.Y     = new double[n+1];
        this.H     = new double[n+1];
        this.A     = new double[n+1];
        this.B     = new double[n+1];
        this.C     = new double[n+1];
        this.Alpha = new double[n+1];
        this.Beta  = new double[n+1];
        this.Solutie = new ArrayList<>();
        X = pct;
        this.a = a1; this.b = b1;
        this.w = f_deriv(a);
        this.z = f_deriv(b);
        this.m = m1;
    }
 
    public void SetareValori()
    {
        for (int i = 1; i <= n; i++)
        {
            Y[i] = f(X[i]);
        }
        
        for (int i = 1; i <= n - 1; i++)
        {
            H[i] = X[i + 1] - X[i];
        }
        
        for (int i = 1; i <= n - 2; i++)
        {
            A[i + 1] = H[i + 1] / (H[i] + H[i + 1]);
            B[i + 1] = H[i] / (H[i] + H[i + 1]);
            C[i + 1] = (3.0 * H[i + 1] * (Y[i + 1] - Y[i]) / (H[i] * (H[i] + H[i + 1]))) + 
                       (3.0 * H[i] * (Y[i + 2] - Y[i + 1]) / (H[i + 1] * (H[i] + H[i + 1])));
        }
        System.out.println();
        System.out.println();
        double precizie = Math.pow(10, -10);
        InitMatrice M = new InitMatrice(10,precizie,A,B,C,w,z);
        V = M.Jacobi(4);
       
        
        for (int i = 1; i <= n - 1; i++)
        {
            Alpha[i] = (3.0 * (Y[i + 1] - Y[i]) / (H[i] * H[i])) - ((V[i + 1] + 2.0 * V[i]) / H[i]);
            
            Beta[i]  = (-2.0 * (Y[i + 1] - Y[i]) / (H[i] * H[i] * H[i])) + ((V[i + 1] + V[i]) / (H[i] * H[i]));
        }
    }
    
    public void RezolvareSpline()
    {
        SetareValori();
        double p = 0;
        double h =0.0;
        double[] t = new double[m+1];
        
        h = (2 * 3.14)/(m);

        for (int j = 0; j <= m-1; j++)
        {
            t[j] = a + h * j;
            
            for (int i = 1; i <= n - 1; i++)
            {
                if (t[j] >= X[i] && t[j] <= X[i + 1])
                {
                    p = Y[i] + V[i] * (t[j] - X[i]) + Alpha[i] * (t[j] - X[i]) * (t[j] - X[i]) +
                            Beta[i] * (t[j] - X[i]) * (t[j] - X[i]) * (t[j] - X[i]);
                    
                    System.out.println("Valoarea functiei f in punctul "+p+" = " 
                    					+ f(t[j])+ " ,unde t["+j+"] = " + t[j]);
                }
            }
                
            System.out.println("Solutia "+j+" = " + p + "\n");
            
            
            Puncte pt = new Puncte(z,p);
            Solutie.add(pt);
            
        }

    }
}

public class Tema5_ex1 {

    private static Scanner sc;

	public static void main(String[] args)
    {
       
        double[] puncte = new double[11];
     
        for (int i = 1; i <= 10; i++)
        {
            puncte[i] = (2.0 * Math.PI * (i - 1)) / 9;
        }
        
        
        sc = new Scanner(System.in);
        System.out.println("m= ");
        int m = sc.nextInt();
        
        Interpolare_spline is = new Interpolare_spline(puncte, 0.0, 6.5, m,10);
        System.out.println("Valorile functiilor spline cubice de interpolare: ");
        is.RezolvareSpline();
        
    }
    
}
