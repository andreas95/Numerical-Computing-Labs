import java.util.Scanner;

import static java.lang.Math.*;

public class Tema4 {

    public static void main(String[] args) {
        double m;
        contractie();
        gaussiedel();
        UnuB();
        System.out.println("\n\n\n\nm=");
        Scanner S=new Scanner(System.in);
        m=S.nextDouble();
        System.out.println("\n\n\n");
        DoiIa(m);
        System.out.println("\n\n\n");
        DoiIIa(m);
        System.out.println("\n\n\n");
        DoiIIb(m);
    }

    static double derv(double x) {
        double dervx=0.0;
        dervx=21*pow(-8*pow(x,3)+11*x+1,2)*(-24*x*x+11)-10*(-24*x*x+11)-1;
        return dervx;
    }

   static double f(double x) {
        double fx=0.0;
        fx=7*pow(-8*pow(x,3)+11*x+1,3)-10*(-8*pow(x,3)+11*x+1)-x+1;
        return fx;
    }

   static void contractie() {
        double y1,y2,x1=0,x2=0,abs1,abs2,maxim;
        int n=0;
        boolean conditie=true;
        while(conditie) {
            n++;
            y1=(1.0/10 * (7*pow(x1,3) - x2 +1));
            y2=(1.0/11 * (8*pow(x2,3) + x1 -1));
            abs1=y1-x1;
            abs2=y2-x2;
            if (abs1<0)
                abs1=abs1* (-1);
            if (abs2<0)
                abs2=abs2* (-1);
            maxim= max(abs1,abs2);
            x1=y1;
            x2=y2;
            conditie = (maxim>=0.0000001);
        }
        System.out.println();
        System.out.println("Metoda Contractiei: "+x1+" "+x2);
    }

    static void gaussiedel()
    {
        double y1,y2,x1=0,x2=0,abs1,abs2,maxim;
        int n=0;
        boolean conditie=true;
        while(conditie) {
            n++;
            y1=(1.0/10 * (7*pow(x1,3) - x2 +1));
            y2=(1.0/11 * (8*pow(x2,3) + y1 -1));

            abs1=y1-x1;
            abs2=y2-x2;
            if (abs1<0)
                abs1=abs1* (-1);
            if (abs2<0)
                abs2=abs2* (-1);
            maxim= max(abs1,abs2);
            x1=y1;
            x2=y2;
            conditie = (maxim>=0.0000001);
        }
        System.out.println();
        System.out.println("Metoda GaussSiedelL "+x1+" "+x2);
    }

    static double rap(double x) {
        return f(x)/derv(x);
    }

    static boolean check(double start , double end, double x) {
        x=x-rap(x);
        if(x<start||x>end)return false;
        return true;
    }

    static double dublu(double x) {
        return -1008*x*pow(-8*pow(x,3)+11*x+1,2)+42*pow(11-24*x*x,2)*(-8*pow(x,3)+11*x+1)+480*x;
    }

    static void newton(double start,double end) {
        boolean test;
        double i,x,y;
        for(i=start;i<end;i=i+0.0001) {
            if(f(i)*f(i+0.0001)<0) {
                x=i;
                while(f(x)*dublu(x)<=0)x=x+0.0001;
                test=check(i,i+0.005,x);
                if(test) {
                    y=x;
                    x=-8*pow(y,3)+11*y+1;
                    System.out.println(x+" "+y);
                }
            }
        }
    }

    static void UnuB() {
        System.out.println();
        System.out.println("Metoda Newton:\n\n\tX\tY\n");
        newton(-2,2);
        System.out.println();
        System.out.println("Metoda Newton simplificata:\n\n\tX\tY\n");
        newton(-2,2);
    }

    static double f1(double x) {
        return 1.0/(1+100*x*x);
    }

    static double interpolare1(double x) {
        double p[]={-1,-0.6666666,-0.2,-0.1,0,0.1,0.2,0.6666666,1};
        int i,j;
        double lx=0,t;
        for(i=0;i<9;i++)
        {
            t=f1(p[i]);//f(x) sau y
            for(j=0;j<i;j++)
                t=t*(x-p[j])/(p[i]-p[j]);
            for(j=i+1;j<9;j++)
                t=t*(x-p[j])/(p[i]-p[j]);
            lx=lx+t;
        }
        return lx;
    }

    static void DoiIa(double m)
    {
        double i;
        System.out.println("Exercitiul 2 punctul 1(A)\n\n");
        for(i=-1.5;i<=1.5;i=i+3.0/(m))
            System.out.print(interpolare1(i)+", ");
        System.out.println();
        System.out.println("Exercitiul 2 punctul 1(B)\n\n");
        for(i=-1.5;i<=1.5;i=i+3.0/(m))
            System.out.print(interpolare1(i)+", ");
    }

    static double f2(double x) {
        return pow(x,10);
    }

    static double interpolare2(double x) {
        double p[]={1,2,3,4,5,6,7,8,9,10};
        int i,j;
        double lx=0,t;
        for(i=0;i<=9;i++)
        {
            t=f2(p[i]);
            for(j=0;j<i;j++)
                t=t*(x-p[j])/(p[i]-p[j]);
            for(j=i+1;j<10;j++)
                t=t*(x-p[j])/(p[i]-p[j]);
            lx=lx+t;
        }
        return lx;
    }

    static void DoiIIa(double m) {
        double i;
        System.out.println("Exercitiul 2 punctul 2(A)\n\n");
        for(i=1;i<10;i=i+9.0/(m-1))
            System.out.print(interpolare2(i)+", ");
        System.out.println();interpolare2(10.0);
    }

    static double interpolare3(double x) {
        double p[]={-10,-9,-8,-7,-6,-5,-4,-3,-2,-1};
        int i,j;
        double lx=0,t;
        for(i=0;i<9;i++)
        {
            t=f2(p[i]);
            for(j=0;j<i;j++)
                t=t*(x-p[j])/(p[i]-p[j]);
            for(j=i+1;j<9;j++)
                t=t*(x-p[j])/(p[i]-p[j]);
            lx=lx+t;
        }
        return lx;
    }

    static void DoiIIb(double m) {
        double i;
        System.out.println("Exercitiul 2 punctul 2(B)\n\n");
        for(i=-10;i<-1;i=i+9.0/(m-1))
            System.out.print(interpolare3(i)+", ");
        System.out.println(interpolare3(-1.0));
    }

}
