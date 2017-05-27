/**
 * @author Ciort Elena Cristina
 */
import java.util.Scanner;

public class Tema5_ex3
{
	public static double X;
	private static Scanner sc;
	static double f(double y)
	{
		return 3.0*Math.pow(X,-3)*(Math.pow(y,3)/(Math.pow(Math.E,y)-1));
	}

static public void formulaDreptunghiului(int n)
{
	double rezultat=0.0;
	double pasi=1.0;
	double a=0.0;
	double b=X;
	
	do
	{
		double h=(b-a)/n;
		double suma=0.0;
		for(int i=1;i<=n;i++)
		{
			suma+=f(a+((2.0*i-1.0)*h)/2.0);
		}
		double sigma=(h)*suma; 
		rezultat=sigma;
		pasi++;
	
	}while(n>pasi);
	System.out.println("Pentru formula dreptunghiului cu "+n+" pasi rezultat="+rezultat);
	System.out.println();
}

static public void formulaNewtonCotes2(int n)
{
	double rezultat=0.0;
	double k=1.0;
	double a=0.0;
	double b=X;
	
	do
	{
		double h=(b-a)/n;
		double suma=0.0;
		for(int i=1;i<=n;i++)
		{
			suma+=f(a+((3.0*i-2.0)*h)/3.0)+f(a+((3.0*i-1.0)*h)/3.0);
		}
		double sigma=(h/2.0)*suma; 
		rezultat=sigma;
		k++;
	
	}while(n>k);
	System.out.println("Pentru formula Newton Cotes pentru n = 2 cu "+n+" pasi rezultat="+rezultat);
	System.out.println();
}

static public void formulaNewtonCotes3(int n)
{
	double rezultat=0.0;
	double pasi=1.0;
	double a=0.0;
	double b=X;
	
	do
	{
		double h=(b-a)/n;
		double suma=0.0;
		for(int i=1;i<=n;i++)
		{
			suma+=2*f(a+((4.0*i-3.0)*h)/4.0)-f(a+((2.0*i-1.0)*h)/2.0)+2*f(a+((4.0*i-1.0)*h)/4.0);
		}
		double sigma=(h/3.0)*suma; 
		rezultat=sigma;
		pasi++;
	
	}while(n>pasi);
	System.out.println("Pentru formula Newton Cotes pentru n = 3 cu "+n+" pasi rezultat="+rezultat);
	System.out.println();
}
static public void formulaNewtonCotes4(int n)
{
	double rezultat=0.0;
	double pasi=1.0;
	double a=0.0;
	double b=X;

	do
	{
		double h=(b-a)/n;
		double suma=0.0;
		for(int i=1;i<=n;i++)
		{
			suma+=11*f(a+((5.0*i-4.0)*h)/5.0)+f(a+((5.0*i-3.0)*h)/5.0)+
					f(a+((5.0*i-2.0)*h)/5.0)+11*f(a+((5.0*i-1.0)*h)/5.0);
		}
		double sigma=(h/24.0)*suma; 
		rezultat=sigma;
		pasi++;
	
	}while(n>pasi);
	System.out.println("Pentru formula Newton Cotes pentru n = 4 cu "+n+" pasi rezultat="+rezultat);
	System.out.println();
}


	public static void main(String[] args)
	{
		Tema5_ex3.X=1.0;
		
		sc = new Scanner(System.in);
		System.out.print("n = ");
		int n = sc.nextInt();
		System.out.println();
		
		formulaDreptunghiului(n);
		formulaNewtonCotes2(n);
		formulaNewtonCotes3(n);
		formulaNewtonCotes4(n);
		
	
	}
}

