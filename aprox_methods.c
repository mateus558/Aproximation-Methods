#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define max(a,b) ((a) > (b) ? (a) : (b))
#define EPS 1E-5
#define MAXIT 1E5

int iterations = 0;

double f(double x){
	return sin(x);
}

double df(double x){
	return cos(x);
}

double itr_f(double x){
	return pow(1 + x, 1.0/3.0);
}

int value_error(double ant, double actual){
	return fabs(f(actual)) > EPS;
}

double val_error(double ant, double actual){
	return fabs(f(actual));
}

int absolute_error(double ant, double actual){
	return fabs(f(actual) - f(ant)) > EPS;
}

double abs_error(double ant, double actual){
	return fabs(actual - ant);
}

double rel_error(double ant, double actual){
	return fabs(actual - ant)/fabs(actual);
}

int relative_error(double ant, double actual){	//k and k+1
	return (fabs(actual - ant) > EPS * max(1, fabs(actual)));
}

double regulafalsi_method(double x0, double x1, int (*stop_criterion)(double, double), double (*error)(double, double)){
	double xk = 2*EPS, x = x0;
	
	if(f(x0)*f(x1) < 0){
		printf("------------------------------------------------------------------------------------\n");
		printf("itr\t   x0\t           x1\t           f(x0)\t          f(x1)\t          m\t          error\n");
		printf("------------------------------------------------------------------------------------\n");
		printf("%i\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", iterations, x0, x1, f(x0), f(x1), x, error(x0, x1));
	
		while(stop_criterion(xk, x) && iterations < MAXIT){
			if(f(x0)*f(x1) < 0){
				xk = x0;
				x0 = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0));
				x = x0;
			}else{
				xk = x1;
				x1 = x1 - f(x1)*(x1 - x0)/(f(x1) - f(x0));	
				x = x1;
			}
			iterations++;
			printf("%i\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", iterations, x0, x1, f(x0), f(x1), x, error(xk, x));
		}
		
		if(iterations >= MAXIT){
			printf("The method did not reached the convergency.\n");
			exit(1);
		}
	}else{
		printf("Error: Invalid interval\n");
		exit(1);
	}
	
	return x;
} 

double quasinewton_method(double x0, double h, int (*stop_criterion)(double, double), double (*error)(double, double)){
	double xk = 2*EPS, d;
	
	printf("---------------------------------------------------------------------\n");
	printf("itr\t   x0\t           f(x)\t          f'(x + h)\t          error\n");
	printf("---------------------------------------------------------------------\n");
	printf("%i\t%lf\t%lf\t%lf\t%lf\n", iterations, x0, f(x0), d, error(xk, x0));
	
	while(stop_criterion(xk, x0) && iterations < MAXIT){
		xk = x0;
		d = (f(x0 + h) - f(x0))/h;
		x0 = x0 - f(x0)/d; 
		iterations++;
		printf("%i\t%lf\t%lf\t%lf\t%lf\n", iterations, x0, f(x0), d, error(xk, x0));
	}
	
	if(iterations >= MAXIT){
		printf("The method did not reached the convergency.\n");
		exit(1);
	}
	
	return xk;
}

double fixedpoint_method(double x0, int (*stop_criterion)(double, double), double (*error)(double, double)){
	double xk = 2*EPS;
	
	printf("--------------------------------------------------------------------\n");
	printf("itr\t   x\t           f(x)\t          itr_f(x)\t          error\n");
	printf("--------------------------------------------------------------------\n");
	printf("%i\t%lf\t%lf\t%lf\t%lf\n", iterations, x0, f(x0), itr_f(x0), error(xk, x0));
	
	while(stop_criterion(xk, x0) && iterations < MAXIT){
		xk = x0;
		x0 = itr_f(x0);
		iterations++;
		printf("%i\t%lf\t%lf\t%lf\t%lf\n", iterations, x0, f(x0), itr_f(x0), error(xk, x0));
	}
	
	if(iterations >= MAXIT){
		printf("The method did not reached the convergency.\n");
		exit(1);
	}
	
	return x0;
}

double secant_method(double x0, double x1, int (*stop_criterion)(double, double), double (*error)(double, double)){
	double xk = x0, nxk = 2*EPS, d;
	
	printf("------------------------------------------------------------------------------------\n");
	printf("itr\t   x0\t           x1\t           f(x)\t          f'(x)\t          error\n");
	printf("------------------------------------------------------------------------------------\n");
	printf("%i\t%lf\t%lf\t%lf\t%lf\t%lf\n", iterations, xk, x1, f(xk), d, error(nxk, xk));
	
	while(stop_criterion(x1, xk) && iterations < MAXIT){
		nxk = xk;
		d = (f(xk) - f(x1))/(xk - x1);
		x1 = xk;
		xk = xk - f(xk)/d; 
		iterations++;
		printf("%i\t%lf\t%lf\t%lf\t%lf\t%lf\n", iterations, xk, x1, f(xk), d, error(nxk, xk));
	}
	
	if(iterations >= MAXIT){
		printf("The method did not reached the convergency.\n");
		exit(1);
	}
	
	return xk;
}

double newton_method(double x0, int (*stop_criterion)(double, double), double (*error)(double, double)){
	double xk = x0, nxk = 2*EPS;
	
	printf("--------------------------------------------------------------------\n");
	printf("itr\t   x\t           f(x)\t          f'(x)\t          error\n");
	printf("--------------------------------------------------------------------\n");
	printf("%i\t%lf\t%lf\t%lf\t%lf\n", iterations, xk, f(xk), df(xk), error(nxk, xk));
	
	while(stop_criterion(nxk, xk) && iterations < MAXIT){
		nxk = xk;
		xk = xk - f(xk)/df(xk); 
		iterations++;
		printf("%i\t%lf\t%lf\t%lf\t%lf\n", iterations, xk, f(xk), df(xk), error(nxk, xk));
	}
	
	if(iterations >= MAXIT){
		printf("The method did not reached the convergency.\n");
		exit(1);
	}
	
	return xk;
}

double bissection_method(double a, double b, int (*stop_criterion)(double, double), double (*error)(double, double)){
	double m = b, ant = 2*EPS;
	
	if(f(a) * f(b) < 0){
	
		printf("-----------------------------------------------------------------------------------------------\n");
		printf("itr\t   x0\t           x1\t           f(x0)\t          f(x1)\t          m\t          error\n");
		printf("-----------------------------------------------------------------------------------------------\n");
		printf("%i\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", iterations, a, b, f(a), f(b), m, error(a, b));
		
		while(stop_criterion(ant, m) && iterations < MAXIT){
			m = (a + b)/2;
			
			if(f(a)*f(m) < 0){
				ant = b;
				b = m;
			}else{
				ant = a;
				a = m;
			}
			
			iterations++;
			printf("%i\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", iterations, a, b, f(a), f(b), m, error(a, b));

		}
		
		if(iterations >= MAXIT){
			printf("The method did not reached the convergency.\n");
			exit(1);
		}
	}else{
		printf("Error: Invalid interval\n");
		exit(1);
	}
	
	return m;
}

int main(){
	double a, b, x, x0;
	int method;

	printf("Select method [0-bissection 1-newton 2-secant 3-fixed_point 4-regula_falsi 5-quasi_newton]: ");	
	scanf("%i",&method);
	
	switch(method){
		case 0:
			printf("a = ");
			scanf("%lf", &a);
			printf("b = ");
			scanf("%lf", &b);
	
			x = bissection_method(a, b, value_error, val_error);
			break;
		case 1:
			printf("x0 = ");
			scanf("%lf", &x0);

			x = newton_method(x0, relative_error, rel_error);
			break;
		case 2:
			printf("x0 = ");
			scanf("%lf", &a);
			printf("x1 = ");
			scanf("%lf", &b);
	
			x = secant_method(a, b, value_error, val_error);
			break;
		case 3:
			printf("x0 = ");
			scanf("%lf", &x0);
			
			x = fixedpoint_method(x0, relative_error, rel_error);
			break;
		case 4:
			printf("x0 = ");
			scanf("%lf", &a);
			printf("x1 = ");
			scanf("%lf", &b);
	
			x = regulafalsi_method(a, b, absolute_error, abs_error);
			break;
		case 5:
			printf("x0 = ");
			scanf("%lf", &a);
			printf("h = ");
			scanf("%lf", &b);
	
			x = quasinewton_method(a, b, relative_error, rel_error);
			break;
		default:
			printf("Invalid method!");
			return 1;
			break;
	}

	printf("\nx = ");
	printf("%f\n", x);
	printf("Number of iterations: ");
	printf("%i\n", iterations+1);
	
	return 0;
}
