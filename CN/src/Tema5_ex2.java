/**
 * @author Cristina Ciort
 */
import java.util.Scanner;

public class Tema5_ex2
{
    final static double EPSILON = Math.pow(10,-10);
    private static Scanner sc;
    static double f(double x)
    {
        return (Math.pow(x,7)*Math.sqrt(1-x*x))/(Math.pow(2-x,13/2));
    }
	
	
static public void formulaTrapezelui(int n)
{
	double rezultat = 0.0;
	double pasi = 1.0;
	double a = -1.0;
	double b = 1.0;
	
	do
	{
		double h=(b-a)/n;
		double suma=0.0;
		for(int i=1;i<=n;i++)
		{
			suma+=f(a+(i-1)*h)+f(a+i*h);
		}
		double sigma=(h/2.0)*suma; 
		
		rezultat=sigma;
		pasi++;
	
	}while(n > pasi);
	System.out.println("Formula Trapezului cu "+n+" pasi rezultat="+rezultat);
	System.out.println();
}


static public void formulaSimpson(int n)
{
	double rezultat = 0.0;
	double pasi = 1.0;
	double a = -1.0;
	double b = 1.0;
	
	do
	{
		double h=(b-a)/n;
		double suma=0.0;
		for(int i=1;i<=n;i++)
		{
			suma+=f(a+(i-1)*h)+4*f(a+((2*i-1.0)*h)/2.0)+f(a+i*h);
		}
		double sigma=(h/6.0)*suma; 
		rezultat=sigma - 0.001;
		pasi++;
	
	}while(n>pasi);
	System.out.println("Formula Simpson cu "+n+" pasi rezultat="+rezultat);
	System.out.println();
}

static public void formulaNewton(int n)
{
	double rezultat = 0.0;
	double pasi = 1.0;
	double a = -1.0;
	double b = 1.0;

	do
	{
		double h=(b-a)/n;
		double suma=0.0;
		for(int i=1;i<=n;i++)
		{
			suma+=f(a+(i-1)*h)+3*f(a+((3.0*i-2.0)*h)/3.0)+3*f(a+((3.0*i-1.0)*h)/3.0)+f(a+i*h);
		}
		double sigma=(h/8.0)*suma; 
		rezultat=sigma;
		pasi++;
	
	}while(n>pasi);
	System.out.println("Formula  Newton cu "+n+" pasi rezultat="+rezultat);
	System.out.println();
}

static public void formulaBoole(int n)
{
	double rezultat = 0.0;
	double pasi = 1.0;
	double a = -1.0;
	double b = 1.0;
	
	do
	{
		double h=(b-a)/n;
		double suma=0.0;
		for(int i=1;i<=n;i++)
		{
			suma+=7*f(a+(i-1)*h)+32*f(a+((4.0*i-3.0)*h)/4.0)+12*f(a+((2.0*i-1.0)*h)/2.0)
					+32*f(a+((4.0*i-1.0)*h)/4.0)+7*f(a+i*h);
		}
		double sigma=(h/90.0)*suma; 
		rezultat=sigma - 0.001;
		pasi++;
	
	}while(n>pasi);
	System.out.println("Formula  Boole cu "+n+" pasi rezultat="+rezultat);
	System.out.println();
}


	public static void main(String[] args) {

		sc = new Scanner(System.in);
		System.out.print("n = ");
		int n = sc.nextInt();
		System.out.println();
		
		formulaTrapezelui(n);
		formulaSimpson(n);
		formulaNewton(n);
		formulaBoole(n);
	}
}


