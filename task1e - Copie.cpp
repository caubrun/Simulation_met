// task1.cpp
//

/*Task 1: Sum of Random Variables. For Sn = L1 + ... + Ln, assume Li are
independent Bernoulli variables with probabilities P(Li = 1) = pi = 1-i/n+1. Use
simulation to generate Sn and compute its probabilities P(Sn = k) and compare the
results with those obtained using Andersen-Sidenius-Basu algorithm and Hull-White
algorithm. For Sn = X1 +...+ Xn, assume Xi are iid random variables (consider
two cases: Bernoulli X1  B(p) and Poisson X1  Pos()). Use simulation to
compute the tail probabilities P(Sn  nx) for x  E[X1] and compare the results
with those from Cramer's theorem and the Central Limit Theorem as n increases.*/

#include <iostream>
#include <random>
#include <vector>
#include <cmath>

using namespace std;

double Andersen_sidenius_Basu(int k, int j, int n) {
	//double P_S_j_K = 0, P_S_j_k;
	if ((k > n) | (j < 0) | (k < 0)) {
		//cout << "0 case : n=" << n << " and k = " << k << endl;
		return 0;
	}
	else {
		if ((j == 0) & (k == 0)) {
			//cout << "1 case" << endl;
			return 1;
		}
		else {
			if (k == j) {
				double Pk = 1;
				for (unsigned int i = 0; i <= j; i++) {
					double P_L = 1 - ((double)i / ((double)n + 1));
					Pk *= P_L;
				}
				return Pk;
			}
			else {
				double P_L = 1 - ((double)j / ((double)n + 1));
				//cout << P_L << "k = " << k << " ; n =" << n << endl;
				return Andersen_sidenius_Basu(k, j - 1, n) * (1 - P_L) + Andersen_sidenius_Basu(k - 1, j - 1, n) * (P_L);
			}
		}
	}
}


double Hull_White(int n, int k) {
	double pi_n = 1;
	for (unsigned int i = 1; i <= n; i++) {
		double PL = 1 - ((double)i / ((double)n + 1));
		pi_n *= (1 - PL);
	}
	if (k == 0) {
		return pi_n;
	}

	double cj;
	vector<double> Vn;
	vector<double> Un;
	for (unsigned int i = 0; i < k; i++) {
		//cout << "k = " << i << endl;
		double Vi = 0;
		for (unsigned int j = 1; j <= n; j++) {
			double PL = 1 - ((double)j / ((double)n + 1));
			cj = PL / (1 - PL);
			//cout << "P = " << PL << endl;
			//cout << "ci = " << ci << endl;
			Vi += pow(cj, i + 1);
		}
		//cout << "i = " << i << "; Vi =" << Vi << endl;
		Vn.push_back(Vi);
		double Ui = 0;
		if (Un.size() == 0) {
			Un.push_back(1);
		}

		for (unsigned int l = 0; l < Vn.size(); l++) {
			//cout << "V" << l << " = " << Vn[l]<<endl;
			//cout << "U" << Un.size() - l-1 << " = " << Un[Un.size() - l-1] << endl;
			Ui += pow(-1, l) * (Vn[l] * Un[Un.size() - l - 1]);
		}
		Ui = Ui / (Un.size());
		Un.push_back(Ui);
	}

	/*for (unsigned int i = 1; i <= n; i++) {
		double PL = 1 - ((double)i / ((double)n + 1));
		pi_n *= (1 - PL);
	}*/
	return pi_n * Un[Un.size() - 1];
}




default_random_engine generator;
uniform_real_distribution<double> distribution(0.0, 1.0);
int generate_Sn(int n) {
	// With C++ tool	
	//double U_i;
	double p_i;
	double L_i;
	double S_n = 0;

	for (unsigned int i = 1; i < n + 1; i++) {
		double U_i = distribution(generator);
		//cout << U_i << endl;
		p_i = 1 - ((double)i / ((double)n + 1));
		//cout << "p : " << p_i << " ; U : " << U_i << endl;
		if (p_i > U_i) {
			L_i = 1;
		}
		else {
			L_i = 0;
		}
		S_n += L_i;
	}
	return S_n;
}

double PSK(int N, int n, int k) {
	vector<int> Sn_values;

	for (unsigned int i = 0; i < N; i++) {
		int S_n = generate_Sn(n);
		//cout << S_n << endl;
		Sn_values.push_back(S_n);
	}
	int count_k = 0;
	for (unsigned int i = 0; i < Sn_values.size(); i++) {
		if (Sn_values[i] == k) {
			count_k += 1;
		}
	}
	double PSk = (double)count_k / (double)N;
	return PSk;
}


default_random_engine generator_poisson;
//poisson_distribution<int> distribution(4.1);


int poisson(double lambda, int n) {
	poisson_distribution<int> poisson_distribution(lambda);
	int S_n = 0;
	for (unsigned int i = 0; i < n; i++) {
		int Xi_poi = poisson_distribution(generator_poisson);
		//cout << Xi_poi << endl;
		S_n += Xi_poi;
	}
	return S_n;
}

double poisson_dist(double x, int N, int n, double lambda) {
	vector<int> Sn_values;
	for (unsigned int i = 0; i < N; i++) {
		int S_n = poisson(lambda, n);
		Sn_values.push_back(S_n);
	}
	int count = 0;
	for (unsigned int k = 0; k < Sn_values.size(); k++) {
		if (Sn_values[k] >= n * x) {
			count += 1;
		}
	}
	return (double)count / (double)N;
}


double cramer(double x, int N, int n, double lambda) {
	return log(poisson_dist(x, N, n, lambda)) / (double)n;
}

double gamma_star(double x, double lambda) {
	return log(x / lambda) * x - x + lambda;

}

double central_limit_poisson(double x, double N, int n, double lambda) {
	double X = (sqrt(lambda) * x ) / sqrt((double)n)+lambda;
	cout << X << endl;
	return poisson_dist(X, N, n, lambda);
}


#include <iomanip>
double normalCDF(double x)
{
	return erfc(-x / sqrt(2)) / 2;
}



int main()
{
	///// Generate n independent variables Ui uniform on [0,1]
	// QQQ : Should we use congrential generator ? C++ tool ? 

	// Distribution of Sn
	int N = 1000; // number of Simulation
	int n = 100; // n
	//for (unsigned int k = 0; k <= n; k++) {
	//	cout << "Result for Simulation : " << PSK(N, n, k) << endl;

	//	//  Andersen-Sidenius-Basu algorithm
	//	// Test
	//	cout << "Result for Andersen Sidenius Basu : " << Andersen_sidenius_Basu(k, n, n) << endl;


	//	// Hull-White algorithm
	//	cout << "Result for Hull White : " << Hull_White(n, k) << endl;
	//}

	double x = 5.02;
	double lambda = 5;
	for (int i = 1; i < 10; i++) {
		n = i * 500;
		cout << poisson_dist(x, N, n, lambda) << endl;
		cout << cramer(x, N, n, lambda) << endl;
		cout << -gamma_star(x, lambda) << endl;
		cout << central_limit_poisson(x, N, n, lambda) << endl;
		cout << 1 - normalCDF(x) << endl;
	}




	return 0;

}

