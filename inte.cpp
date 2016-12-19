#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <cmath>

using namespace std;

typedef vector<vector<double> > dMatrix;
typedef vector<vector<int> > iMatrix;
struct Point{
	double x, y;
};

bool comp(Point a, Point b);
double f(double x);
double g(double x);
double f1(double x);
double g1(double x);
double lagrange_interpol(vector<Point> points, double x);
bool load_points(string fname, vector<Point> &points);
dMatrix swap_lines(dMatrix A, int N, int L1, int L2);
vector<double> gauss_elimination(dMatrix A, int N);
vector<double> LLSR_generalized(vector<Point> points, int degree, double (*p[]) (double x),double &quad_mean);
vector<double> LLSR_polynomial(vector<Point> points, int degree, double &quad_mean);
vector<double> nLLSR_exponential(vector<Point> points, double &quad_mean);
vector<double> nLLSR_four(vector<Point> points, double &quad_mean);
 
double (*q[2]) (double x) = {
	f, g
}; 

double (*p[2]) (double x) = {
	f1, g1
};

int main(){
	int i, d, o;
	double x, quad_mean = 0.0;
	vector<double> c;
	vector<Point> points;
	
	load_points("in", points);
	
	sort(points.begin(), points.end(), comp);
	
	cout << "Loaded points: ";
	for(Point p: points)
		cout << "(" << p.x << ", " << p.y << ") | ";
	cout << endl;
	
	cout << "\n0 - Evaluate a point by lagrange interpolation; 1..5 - Least squares method functions" << endl;  
	cout << "Select an option: ";
	cin >> o;
	
	switch(o){
		case 0:
			cout << "Point to evaluate: ";
			cin >> x;
			
			if(x < points[0].x || x > points[points.size()-1].x){
				cout << "Point out of interpolation range." << endl;
				return 0;
			}
			
			x = lagrange_interpol(points, x);
			cout << "Value: " << x << endl;
			return 0;
			break;
		case 1:
			cout << "\nLeast square method with n-degree polynomial\n" << endl;
			cout << "Degree of the polynomial: ";
			cin >> d;
			cout << endl;
			
			c = LLSR_polynomial(points, d, quad_mean);
			break;
		case 2:
			cout << "\nLeast square method with exponential\n" << endl;
			
			c = nLLSR_exponential(points, quad_mean);		
			break;
		case 3:
			cout << "\nLeast square method with exponential\n" << endl;
			
			c = nLLSR_four(points, quad_mean);		
			break;
		case 4:
			cout << "\nGeneralized Least Square Method 1\n" << endl;
			cout << "f(x) = e^x and g(x) = e^-x\n" << endl;
			
			c = LLSR_generalized(points, 1, q, quad_mean);		
			break;
		case 5:
			cout << "\nGeneralized Least Square Method 2\n" << endl;
			cout << "f(x) = 1/(1 + e^-x) and g(x) = e^-(x^2)\n" << endl;
			
			c = LLSR_generalized(points, 1, p, quad_mean);		
			break;
		default:
			cout << "Unknown option, aborting..." << endl;
			return 0;  
	}
	
	cout << "Quadratic mean error: " << quad_mean << endl;
	cout << endl;
	for(i = 0; i < c.size(); ++i){
		cout << "c[" << i << "]: " << c[i] << endl; 
	}
}

bool comp(Point a, Point b){
	return a.x < b.x;
}

bool load_points(string fname, vector<Point> &points){
	int N, i;
	ifstream input(fname.c_str());
	
	if(!input){
		cout << "File couldn't be opened!" << endl;
		return 0;
	}
	
	input >> N;
	points.resize(N);
	
	for(i = 0; i < N; ++i){
		input >> points[i].x >> points[i].y;
	}
	
	return 1;	
}

double f(double x){
	return exp(x);
}

double g(double x){
	return exp(-x);
}

double f1(double x){
	return 1/(1 + exp(-x));
}

double g1(double x){
	return exp(-x * x);
}

double lagrange_interpol(vector<Point> points, double x){
	int j, k, size = points.size();
	double  y, sum = 0.0, d, sum1 = 0.0;
	vector<double> w(size, 1.0), l(size, 1.0);
	
	/*for(j = 0; j < size; ++j){
		for(k = 0; k < j; ++k){
			w[j] *= (points[j].x - points[k].x);
		}
		for(k = j + 1; k < size; ++k){
			w[j] *= (points[j].x - points[k].x);
		}
		w[j] = 1 / w[j];
	}
	for(j = 0; j < size; ++j){
		l *= (x - points[j].x);
	}
	for(j = 0; j < size; ++j){			
		sum += (w[j] / (x - points[j].x)) * points[j].y;
	}*/
	for(j = 0; j < size; ++j){
		for(k = 0; k < j; ++k){
			l[j] *= (x - points[k].x) / (points[j].x - points[k].x);	
		}
		for(k = j + 1; k < size; ++k){
			l[j] *= (x - points[k].x) / (points[j].x - points[k].x);	
		}	
		sum += points[j].y * l[j];
	}
	
	return sum;
}

dMatrix swap_lines(dMatrix A, int N, int L1, int L2){
	int i, j, k;
	double sum = 0.0, temp; 
	iMatrix P(N, vector<int>(N));
	dMatrix res(N, vector<double>(N));
	vector<double> b(N);
	
	for(i = 0; i < N; ++i){ 
		for(j = 0; j < N; ++j)
			P[i][j] = 0;
		P[i][i] = 1;
	}
	
	for(i = 0; i < N; ++i){
		b[i] = A[i][N];
	}
	
	P[L1][L1] = 0;
	P[L1][L2] = 1;
	P[L2][L2] = 0;
	P[L2][L1] = 1;
	
	for(i = 0; i < N; i++){
    	for(j = 0; j < N; j++){
    		for(k = 0, sum = 0.0; k < N; k++){
    			sum += P[i][k] * A[k][j];
        	}
        	res[i][j] = sum;
    	}
	}
	
	A[L1][N] = b[L2];
	A[L2][N] = b[L1];
	
	for(i = 0; i < N; i++){
    	for(j = 0; j < N; j++){
			A[i][j] = res[i][j];
		}
	}
	
	return A;
}

vector<double> gauss_elimination(dMatrix A, int N){
	int i, j, k, l = -1, flag = 1;
	double ratio, sum = 0.0; 
	vector<double> b(N, 0.0), x(N, 0.0);
	
	for(i = 0; i < N; ++i) b[i] = A[i][N];
	
	for(k = 0; k < N - 1; ++k){
		for(i = k + 1; i < N; ++i){
			if(A[k][k] == 0){
				for(j = 0; j < k; ++j){
					if(l != -1 && A[j][k] != 0 && fabs(A[j][k]) > fabs(A[l][k])) l = j;
				}
				
				for(j = k+1; j < N; ++j){
					if(l != -1 && A[j][k] != 0 && fabs(A[j][k]) > fabs(A[l][k])) l = j;
				}
				
				if(l != -1)
					A = swap_lines(A, N, k, l);
				l = -1;
			}
			ratio = A[i][k]/A[k][k];
			A[i][k] -= A[i][k];
			for(j = k + 1; j < N; ++j){
				A[i][j] -= ratio * A[k][j];
			}
			b[i] -= ratio * b[k];
		}
	}
	
	x[N-1] = b[N-1]/A[N-1][N-1];
	
	for(i = N-1; i >= 0; --i){
		for(j = i+1, sum = b[i]; j < N; ++j){
			sum -= A[i][j]*x[j];
		}
		x[i] = sum/A[i][i];
	}
	
	return x;
}

vector<double> LLSR_polynomial(vector<Point> points, int degree, double &quad_mean){
	int i, j, k, size = degree+1, psize = points.size();
	double quad_rest = 0.0, s = 0.0;
	vector<double> x_sums(2*degree+1, 0.0), y_sums(degree+1, 0.0), c(degree+1);
	dMatrix matrix(size, vector<double>(size + 1, 0.0));
	
	for(i = 0; i < 2*degree+1; ++i){
		for(k = 0; k < psize; ++k){
			x_sums[i] += pow(points[k].x, i); 
			if(i <= degree){
				y_sums[i] += pow(points[k].x, i) * points[k].y;				
			}
		}
	}
	
	for(i = 0; i < size; ++i){
		for(j = 0; j < size; ++j){
			matrix[i][j] = x_sums[i + j];
		}
		matrix[i][j] = y_sums[i];
	}
	
	c = gauss_elimination(matrix, size);
	
	for(i = 0; i < psize; ++i){
		s = points[i].y;
		for(j = size - 1; j >= 0; --j){
			s -= c[j] * pow(points[i].x, j); 
		}
		quad_rest += s * s;
	}
	quad_mean = quad_rest / psize;
	
	return c;
}

vector<double> nLLSR_exponential(vector<Point> points, double &quad_mean){
	int degree = 1;
	int i, j, k, size = degree+1, psize = points.size();
	double quad_rest = 0.0, s = 0.0;
	vector<double> x_sums(2*degree+1, 0.0), y_sums(degree+1, 0.0), c(degree+1);
	dMatrix matrix(size, vector<double>(size + 1, 0.0));
	
	for(i = 0; i < 2*degree+1; ++i){
		for(k = 0; k < psize; ++k){
			x_sums[i] += pow(points[k].x, i); 
			if(i <= degree){
				y_sums[i] += pow(points[k].x, i) * log(points[k].y);				
			}
		}
	}
	
	for(i = 0; i < size; ++i){
		for(j = 0; j < size; ++j){
			matrix[i][j] = x_sums[i + j];
		}
		matrix[i][j] = y_sums[i];
	}
	
	c = gauss_elimination(matrix, size);	
	c[0] = exp(c[0]);
	
	for(i = 0; i < psize; ++i){
		s = log(points[i].y);
		s -= log(fabs(c[1])) + c[0] * points[i].x;
		quad_rest += s * s;
	}
	quad_mean = quad_rest / psize;
	
	return c;
}

vector<double> nLLSR_four(vector<Point> points, double &quad_mean){
	int degree = 1;
	int i, j, k, size = degree+1, psize = points.size();
	double quad_rest = 0.0, s = 0.0;
	vector<double> x_sums(2*degree+1, 0.0), y_sums(degree+1, 0.0), c(degree+1);
	dMatrix matrix(size, vector<double>(size + 1, 0.0));
	
	for(i = 0; i < 2*degree+1; ++i){
		for(k = 0; k < psize; ++k){
			x_sums[i] += pow(points[k].x, i); 
			if(i <= degree){
				y_sums[i] += pow(points[k].x, i) * log(points[k].y) - pow(points[k].x, i) * log(points[k].x);				
			}
		}
	}
	
	for(i = 0; i < size; ++i){
		for(j = 0; j < size; ++j){
			matrix[i][j] = x_sums[i + j];
		}
		matrix[i][j] = y_sums[i];
	}
	
	c = gauss_elimination(matrix, size);	
	
	for(i = 0; i < psize; ++i){
		s = log(points[i].y);
		s -= log(fabs(c[1])) + log(fabs(points[i].x)) + c[0] * points[i].x; 
		quad_rest += s * s;
	}
	quad_mean = quad_rest / psize;
	c[0] = exp(c[0]);
	
	return c;
}

vector<double> LLSR_generalized(vector<Point> points, int degree, double (*p[]) (double x), double &quad_mean){
	int i, j, k, size = degree+1, psize = points.size();
	double quad_rest = 0.0, s = 0.0;
	vector<double> dots(2*degree+1, 0.0), dot_y(degree+1, 0.0), c(degree+1);
	dMatrix u(size, vector<double>(psize));
	dMatrix matrix(size, vector<double>(size + 1, 0.0));
	
	for(i = 0; i < size; ++i){
		for(j = 0; j < psize; ++j){
			u[i][j] = (*p[i]) (points[j].x);
		}
	}
	
	for(i = 0; i < size; ++i){
		for(j = i; j < size; ++j){
			for(k = 0; k < psize; ++k)
				dots[i + j] += u[i][k] * u[j][k];
		}
		for(k = 0; k < psize; ++k){
			dot_y[i] += u[i][k] * points[k].y;
		}
	}

	for(i = 0; i < size; ++i){
		for(j = 0; j < size; ++j){
			matrix[i][j] = dots[i + j];
		}
		matrix[i][j] = dot_y[i];
	}
	
	c = gauss_elimination(matrix, size);
	
	for(i = 0; i < psize; ++i){
		s = points[i].y;
		for(j = size - 1; j >= 0; --j){
			s -= c[j] * u[j][i]; 
		}
		quad_rest += s * s;
	}
	quad_mean = quad_rest / psize;

	return c;
}
