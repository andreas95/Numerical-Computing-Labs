import static java.lang.Math.*;

public class Tema1 {

    private static double err=0.000000001;
    private static double sol;

    public static void main(String[] args) {
      //Metoda bisectiei
        System.out.println("Metoda bisectie:");
        metoda_bisectiei(0,3);
/*
        //Metoda falsei pozitii
        System.out.println("Metoda falsei pozitii:");
        metoda_falsi(0,3);*/

        //Meotada lui Newton
        System.out.println("Metoda lui Newton:");
        metoda_newton();

    /*    //Metoda Secantei
        System.out.println("Metoda Secantei:");
        metoda_secantei(0,0.25);
        metoda_secantei(0.25,0.75);
        metoda_secantei(0.5,1);
        metoda_secantei(0.5,1.35);
        metoda_secantei(0.5,1.65);
        metoda_secantei(0.5,2.85);
        metoda_secantei(0.5,2.15);
        metoda_secantei(0.5,0.95);
        metoda_secantei(0.5,2.65);
        metoda_secantei(0.5,3);

        //Metoda Coardei
        System.out.println("Metoda coardei:");
        metoda_coardei(0,0.25);
        metoda_coardei(0.25,0.75);
        metoda_coardei(0.5,1);
        metoda_coardei(0.5,1.35);
        metoda_coardei(0.5,1.65);
        metoda_coardei(0.5,2.85);
        metoda_coardei(0.5,2.15);
        metoda_coardei(0.5,0.95);
        metoda_coardei(0.5,2.65);
        metoda_coardei(0.5,3);*/

    }

    public static double f(double x) {
        return exp(-1*x)*4*sin(6*x*sqrt(x))-0.2;
    }

    public static double fd(double x) {
        return 36*exp(-x)*sqrt(x)*cos(6*pow(x,(3/2))) - 4*exp(-x)*sin(6*pow(x,(3/2)));
    }

    public static void metoda_bisectiei(double a, double b) {
        double x1,x2,x;
        for (double i=0; i<b; i+=0.25) {
            x1=i;
            x2=i+0.25;
            if (f(x1)*f(x2)<=0) {
                while (x2-x1>err) {
                    x=(x1+x2)/2;
                    if (f(x2)*f(x)<=0) {
                        x1=x;
                    } else {
                        x2=x;
                    }
                }

                if (f(x1)==0) {
                    System.out.printf("%.2f",x1);
                } else if (f(x2)==0) {
                    System.out.printf("%.2f",x2);
                } else
                System.out.printf("%.2f",(x1+x2)/2);
                System.out.println();
            }
        }
    }

    public static void metoda_falsi(double a, double b) {
        double x1, x2, x=0;
        for (double i = 0; i < b; i += 0.25) {
            x1 = i;
            x2 = i + 0.25;
            if (f(x1) * f(x2) < 0) {
                for (int j = 0; j < 100; j++) {
                    x = (x1 * f(x2) - x2 * f(x1)) / (f(x2) - f(x1));
                    if (f(x) == 0) {
                        break;
                    } else if (f(x1) * f(x) < 0) {
                        x2 = x;
                    } else {
                        x1 = x;
                    }
                }
                System.out.printf("%.2f",x);
                System.out.println();
            }
        }
    }

    public static void metoda_newton() {

        System.out.println();
    }

    public static void metoda_secantei(double a,double b) {
        double e=0.0001;
        if (f(a)*f(b)<0) {
            double c;
            do {
                c = (a * f(b) - b * f(a)) / (f(b) - f(a));
                a = b;
                b = c;
            } while (abs(f(c))>0.005);
            System.out.printf("%.2f",c);
            System.out.println();
        }
    }

    public static void metoda_coardei(double a, double b) {
        metoda_secantei(a,b);
        /*double e=0.0001;
        if (f(a)*f(b)<0) {
            double c,x0,x1;
            c=b;
            x0=a;
            do {
                x1 = (c * f(x0) - x0 * f(c)) / (f(x0) - f(c));
                x0=x1;
            } while (abs(f(x1))>0.005);
            System.out.printf("%.2f",x1);
            System.out.println();
        }*/
    }

}
