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

#include <time.h>
#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <iomanip>
#include <functional> 
#include <fstream>
#include<string>
#include <ctime>
#include <ratio>
#include <chrono>
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////
//									  Timer functions 										  //
////////////////////////////////////////////////////////////////////////////////////////////////


// GLOBAL VARIABLE TO STORE STARTING TIME:
//chrono::high_resolution_clock::time_point start_time;
chrono::system_clock::time_point start_time;

/*
void start_timing(chrono::high_resolution_clock::time_point t) {
	start_time = t;
}
*/

double run_time_sec() {
	chrono::system_clock::time_point t2 = chrono::system_clock::now();
	chrono::duration<double> time_span = chrono::duration_cast<std::chrono::microseconds>(t2 - start_time);
	//chrono::duration_cast<chrono::duration<double>>(t2 - start_time);
	double run_duration = time_span.count();
	return run_duration;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//									Global functions		      							  //
////////////////////////////////////////////////////////////////////////////////////////////////
double normalCDF(double x)
{/*normal cdf*/
	return erfc(-x / sqrt(2)) / 2;
}

//Normal PDF
double normalPDF(double x) {
	return (1 / sqrt(2 * M_PI) * exp(-pow(x, 2) * 0.5));
}

////////////////////////////////////////////////////////////////////////////////////////////////
//							  Generation of random numbers      							  //
////////////////////////////////////////////////////////////////////////////////////////////////


class random_gen {
	/*This class is the random generator for the several tasks of our program. In fact, it allows to generate a vecotr of N uniform random variables,
	and to derive N normal variables, N Bernoulli variables, N Poisson variables*/
	unsigned int N;
	vector<double> unif;
public:
	random_gen(int nb_variable) {
		//cout << "Constructing" << endl;
		N = nb_variable;
		int a = pow(7, 5);
		int m = pow(2, 31) - 1;
		time_t t;
		/*int n0 = time(&t);*/
		random_device rd;  //Will be used to obtain a seed for the random number engine
		srand(rd());
		int n0 = rand() % time(&t);
		for (unsigned int i = 0; i < N; i++) {
			n0 = (n0 * a);

			if (n0 < 0) {
				n0 = abs(-m - n0);
			}
			unif.push_back((double)n0 / (double)m);
		}
	}
	void print() {
		for (unsigned int i = 0; i < N; i++) {
			cout << unif[i] << endl;
		}
	}
	vector<double> get_unif() {
		return unif;
	}
	vector<double> normal() {
		//Algorithm (Marsaglia Polar) 
		vector<double> normal_gen;
		while ((int)normal_gen.size() < N) {
			random_gen U(N);

			for (int i = 0; i < N - 1; i++) {
				double v1 = 2 * U.unif[i] - 1;
				double v2 = 2 * U.unif[i + 1] - 1;
				double W = v1 * v1 + v2 * v2;
				if (W < 1) {
					normal_gen.push_back(sqrt(-2 * log(W) / W) * v1);
					normal_gen.push_back(sqrt(-2 * log(W) / W) * v2);
				}
				if (normal_gen.size() == N) {
					break;
				}
			}
		}
		return normal_gen;
	}
	vector<int> poisson(double lambda) {
		/*algorithm poisson random number (Knuth)*/
		double L = exp(-lambda);
		vector<int> poisson_gen;
		for (unsigned int i = 0; i < N; i++) {
			int k = 0;
			double p = 1;
			while (L < p) {
				k += 1;
				// If we build only one uniform variable, it gives always the same, since the clock doesnt change in small time; this is a shortcut
				random_gen U(i +k );
				double u = U.unif[i+k-1];
				p *= u;
			}
			poisson_gen.push_back(k - 1);
		}
		return poisson_gen;
	}
	vector<int> bernoulli(double p) {
		vector<int> bernoulli_gen;
		for (unsigned int i = 0; i < N; i++) {
			if (unif[i] < p) {
				bernoulli_gen.push_back(1);
			}
			else {
				bernoulli_gen.push_back(0);
			}
		}
		return bernoulli_gen;
	}

};


////////////////////////////////////////////////////////////////////////////////////////////////
//						Global task : Call Option Statisitcs								  //
////////////////////////////////////////////////////////////////////////////////////////////////

class Call_EU {
	/*This class allows to compute the global task. In fact, it computes via Monte Carlo Method : the price of the call option,
	the Greeks (likehood ratio and pathwise), the different statistics*/
	double S0, K, r, T, sigma;
	double d1 = (1 / (sigma * sqrt(T)) * (log(S0 / K) + (r + 0.5 * sigma * sigma) * T));
	double d2 = d1 - sigma * sqrt(T);


	vector<double> S_T;
	vector<double> payoff_vec;
	double payoff_mc = 0.0;
	double delta_mc_pw = 0.0;
	double delta_mc_lr = 0.0;
	double vega_mc_pw = 0.0;
	double vega_mc_lr = 0.0;
	double rho_mc_pw = 0.0;
	double rho_mc_lr = 0;
	double gamma_mc_lr_lr = 0.0;
	double gamma_mc_lr_pw = 0.0;
	double gamma_mc_pw_lr = 0.0;
	double sum_square_price = 0.0;
	double sum_square_delta_pw = 0.0;
	double sum_square_vega_pw = 0.0;
	double sum_square_rho_pw = 0.0;
	double sum_square_delta_lr = 0.0;
	double sum_square_vega_lr = 0.0;
	double sum_square_rho_lr = 0.0;
	double sum_square_gamma_lr_lr = 0.0;
	double sum_square_gamma_lr_pw = 0.0;
	double sum_square_gamma_pw_lr = 0.0;
	unsigned long long int m = 0;
	double time_taken = 0;


public:
	Call_EU(double s, double k, double R, double t, double sig) : S0(s), K(k), r(R), T(t), sigma(sig) {};

	void MonteCarlo(unsigned long long int M) {
		S_T.clear();
		payoff_vec.clear();
		random_gen random(M);
		auto start = chrono::high_resolution_clock::now();

		vector<double> normal = random.normal();
		m = M;

		for (int i = 0; i < M; i++) {
			double Z = normal[i];
			double S_T_i = S0 * exp((r - 0.5 * sigma * sigma) * T + sigma * sqrt(T) * Z);

			S_T.push_back(S_T_i);

			if (S_T_i > K) {
				payoff_vec.push_back(S_T_i - K);
				payoff_mc += (S_T_i - K);

				delta_mc_pw += S_T_i / S0;
				delta_mc_lr += ((S_T_i - K) * Z) / (S0 * sigma * sqrt(T));

				vega_mc_pw += (Z * sqrt(T) - sigma * T) * S_T_i;
				vega_mc_lr += (S_T_i - K) * ((Z * Z / sigma) - (Z * sqrt(T)) - (1 / sigma));

				gamma_mc_lr_lr += (S_T_i - K) * ((pow(Z, 2) - 1) / (pow(S0, 2) * pow(sigma, 2) * T) - Z / (pow(S0, 2) * sigma * sqrt(T)));
				gamma_mc_lr_pw += K * Z / (S0 * S0 * sigma * sqrt(T));
				gamma_mc_pw_lr += (S_T_i / (S0 * S0)) * ((Z / (sigma * sqrt(T)) - 1));

				rho_mc_pw += K * T;
				rho_mc_lr += (S_T_i - K) * ((Z * sqrt(T) / sigma) - T);

				sum_square_price += pow((S_T_i - K), 2);
				sum_square_delta_pw += pow((S_T_i / S0), 2);
				sum_square_vega_pw += pow(((Z * sqrt(T) - sigma * T) * S_T_i), 2);
				sum_square_rho_pw += pow(S_T_i * T, 2);

				sum_square_gamma_lr_pw += pow(K * Z / (S0 * S0 * sigma * sqrt(T)), 2);
				sum_square_gamma_lr_lr += pow((S_T_i - K) * ((pow(Z, 2) - 1) / (pow(S0, 2) * pow(sigma, 2) * T) - Z / (pow(S0, 2) * sigma * sqrt(T))), 2);
				sum_square_gamma_pw_lr += pow((S_T_i / (S0 * S0)) * ((Z / (sigma * sqrt(T)) - 1)), 2);

				sum_square_delta_lr += pow(((S_T_i - K) * Z) / (S0 * sigma * sqrt(T)), 2);
				sum_square_vega_lr += pow((S_T_i - K) * ((Z * Z / sigma) - (Z * sqrt(T)) - (1 / sigma)), 2);
				sum_square_rho_lr += pow((S_T_i - K) * ((Z * sqrt(T) / sigma) - T), 2);
			}
			else {
				payoff_vec.push_back(0);
			}
		}

		payoff_mc *= exp(-r * T) / M;
		delta_mc_pw *= exp(-r * T) / M;
		delta_mc_lr *= exp(-r * T) / M;
		vega_mc_pw *= exp(-r * T) / M;
		vega_mc_lr *= exp(-r * T) / M;
		rho_mc_pw *= exp(-r * T) / M;
		rho_mc_lr *= exp(-r * T) / M;
		gamma_mc_lr_lr *= exp(-r * T) / M;
		gamma_mc_lr_pw *= exp(-r * T) / M;
		gamma_mc_pw_lr *= exp(-r * T) / M;
		sum_square_price *= exp(-2 * r * T) / (M * (M - 1));
		sum_square_delta_pw *= exp(-2 * r * T) / (M * (M - 1));
		sum_square_vega_pw *= exp(-2 * r * T) / (M * (M - 1));
		sum_square_rho_pw *= exp(-2 * r * T) / (M * (M - 1));
		sum_square_delta_lr *= exp(-2 * r * T) / (M * (M - 1));
		sum_square_vega_lr *= exp(-2 * r * T) / (M * (M - 1));
		sum_square_rho_lr *= exp(-2 * r * T) / (M * (M - 1));
		sum_square_gamma_lr_lr *= exp(-2 * r * T) / (M * (M - 1));
		sum_square_gamma_lr_pw *= exp(-2 * r * T) / (M * (M - 1));
		sum_square_gamma_pw_lr *= exp(-2 * r * T) / (M * (M - 1));

		auto end = chrono::high_resolution_clock::now();
		time_taken = chrono::duration_cast<chrono::milliseconds>(end - start).count();
	}

	double price_BS() {
		return S0 * normalCDF(d1) - K * exp(-r * T) * normalCDF(d2);
	}

	double delta_bs() {
		return normalCDF(d1);
	}

	double gamma_bs() {
		return (1 / (sigma * S0 * sqrt(T))) * normalPDF(d1);
	}

	double vega_bs() {
		return S0 * sqrt(T) * normalPDF(d1);
	}

	double rho_bs() {
		return K * T * exp(-r * T) * normalCDF(d2);
	}

	double price_mc() {
		return payoff_mc;
	}

	double bias() {
		return price_BS() - price_mc();
	}

	double variance_price() {
		return sum_square_price - (pow(payoff_mc, 2) / (m - 1));
	}

	double variance_delta_pw() {
		return sum_square_delta_pw - (pow(delta_mc_pw, 2) / (m - 1));
	}

	double variance_vega_pw() {
		return sum_square_vega_pw - (pow(vega_mc_pw, 2) / (m - 1));
	}

	double variance_rho_pw() {
		return sum_square_rho_pw - (pow(rho_mc_pw, 2) / (m - 1));
	}

	double variance_delta_lr() {
		return sum_square_delta_lr - (pow(delta_mc_lr, 2) / (m - 1));
	}

	double variance_vega_lr() {
		return sum_square_vega_lr - (pow(vega_mc_lr, 2) / (m - 1));
	}

	double variance_rho_lr() {
		return sum_square_rho_lr - (pow(rho_mc_lr, 2) / (m - 1));
	}

	double variance_gamma_lr_lr() {
		return sum_square_gamma_lr_lr - (pow(gamma_mc_lr_lr, 2) / (m - 1));
	}

	double variance_gamma_lr_pw() {
		return sum_square_gamma_lr_pw - (pow(gamma_mc_lr_pw, 2) / (m - 1));
	}

	double variance_gamma_pw_lr() {
		return sum_square_gamma_pw_lr - (pow(gamma_mc_pw_lr, 2) / (m - 1));
	}

	double get_delta_mc_pw() {
		return delta_mc_pw;
	}

	double get_delta_mc_lr() {
		return delta_mc_lr;
	}

	double get_vega_mc_pw() {
		return vega_mc_pw;
	}

	double get_vega_mc_lr() {
		return vega_mc_lr;
	}

	double get_rho_mc_pw() {
		return rho_mc_pw;
	}

	double get_rho_mc_lr() {
		return  rho_mc_lr;
	}

	double gamma_lr_lr() {
		return gamma_mc_lr_lr;
	}

	double gamma_lr_pw() {
		return gamma_mc_lr_pw;
	}

	double gamma_pw_lr() {
		return gamma_mc_pw_lr;
	}

	double time() {
		return time_taken;
	}

	ostream& Statistics(ostream& os) {

		os << "Number of Simulations: " << m << endl;
		os << "Computation Time: " << time_taken << endl;
		os << "" << endl;

		os << "Statistics of Call Option: Monte Carlo - Price:" << endl;
		os << "Estimate of Call Option Price: " << price_mc() << endl;
		os << "Actual Call Option Price " << price_BS() << endl;
		os << "Bias of estimator: " << bias() << endl;
		os << "Variance of estimator: " << variance_price() << endl;
		os << "Mean Square Error: " << variance_price() + pow(bias(), 2) << endl;
		os << "99.7% Confidence Interval [" << price_mc() - 3 * sqrt(variance_price()) << ", " << price_mc() + 3 * sqrt(variance_price()) << "]" << endl;
		os << "" << endl;

		//Pathwise Statistics
		os << "Statistics of Call Option: Monte Carlo - Delta(PW):" << endl;
		os << "Estimate of Call Option Delta: " << delta_mc_pw << endl;
		os << "Actual Call Option Delta " << delta_bs() << endl;
		os << "Bias of estimator: " << delta_bs() - delta_mc_pw << endl;
		os << "Variance of estimator: " << variance_delta_pw() << endl;
		os << "Mean Square Error: " << variance_price() + pow((delta_bs() - delta_mc_pw), 2) << endl;
		os << "99.7% Confidence Interval [" << delta_mc_pw - 3 * sqrt(variance_delta_pw()) << ", " << delta_mc_pw + 3 * sqrt(variance_delta_pw()) << "]" << endl;
		os << "" << endl;

		os << "Statistics of Call Option: Monte Carlo - Delta(LR):" << endl;
		os << "Estimate of Call Option Delta: " << delta_mc_lr << endl;
		os << "Actual Call Option Delta " << delta_bs() << endl;
		os << "Bias of estimator: " << delta_bs() - delta_mc_lr << endl;
		os << "Variance of estimator: " << variance_delta_lr() << endl;
		os << "Mean Square Error: " << variance_price() + pow((delta_bs() - delta_mc_lr), 2) << endl;
		os << "99.7% Confidence Interval [" << delta_mc_lr - 3 * sqrt(variance_delta_lr()) << ", " << delta_mc_lr + 3 * sqrt(variance_delta_lr()) << "]" << endl;
		os << "" << endl;

		os << "Statistics of Call Option: Monte Carlo - Vega(PW):" << endl;
		os << "Estimate of Call Option Vega: " << vega_mc_pw << endl;
		os << "Actual Call Option Vega " << vega_bs() << endl;
		os << "Bias of estimator: " << vega_bs() - vega_mc_pw << endl;
		os << "Variance of estimator: " << variance_vega_pw() << endl;
		os << "Mean Square Error: " << variance_vega_pw() + pow((vega_bs() - vega_mc_pw), 2) << endl;
		os << "99.7% Confidence Interval [" << vega_mc_pw - 3 * sqrt(variance_vega_pw()) << ", " << vega_mc_pw + 3 * sqrt(variance_vega_pw()) << "]" << endl;
		os << "" << endl;

		os << "Statistics of Call Option: Monte Carlo - Vega(LR):" << endl;
		os << "Estimate of Call Option Vega: " << vega_mc_lr << endl;
		os << "Actual Call Option Vega " << vega_bs() << endl;
		os << "Bias of estimator: " << vega_bs() - vega_mc_lr << endl;
		os << "Variance of estimator: " << variance_vega_lr() << endl;
		os << "Mean Square Error: " << variance_vega_lr() + pow((vega_bs() - vega_mc_lr), 2) << endl;
		os << "99.7% Confidence Interval [" << vega_mc_lr - 3 * sqrt(variance_vega_lr()) << ", " << vega_mc_lr + 3 * sqrt(variance_vega_lr()) << "]" << endl;
		os << "" << endl;

		os << "Statistics of Call Option: Monte Carlo - Rho(PW):" << endl;
		os << "Estimate of Call Option Rho: " << get_rho_mc_pw() << endl;
		os << "Actual Call Option Rho " << rho_bs() << endl;
		os << "Bias of estimator: " << rho_bs() - get_rho_mc_pw() << endl;
		os << "Variance of estimator: " << variance_rho_pw() << endl;
		os << "Mean Square Error: " << variance_rho_pw() + pow((rho_bs() - get_rho_mc_pw()), 2) << endl;
		os << "99.7% Confidence Interval [" << get_rho_mc_pw() - 3 * sqrt(variance_rho_pw()) << ", " << get_rho_mc_pw() + 3 * sqrt(variance_rho_pw()) << "]" << endl;
		os << "" << endl;

		os << "Statistics of Call Option: Monte Carlo - Rho(LR):" << endl;
		os << "Estimate of Call Option Rho: " << get_rho_mc_lr() << endl;
		os << "Actual Call Option Rho " << rho_bs() << endl;
		os << "Bias of estimator: " << rho_bs() - get_rho_mc_lr() << endl;
		os << "Variance of estimator: " << variance_rho_lr() << endl;
		os << "Mean Square Error: " << variance_rho_lr() + pow((rho_bs() - get_rho_mc_lr()), 2) << endl;
		os << "99.7% Confidence Interval [" << get_rho_mc_lr() - 3 * sqrt(variance_rho_lr()) << ", " << get_rho_mc_lr() + 3 * sqrt(variance_rho_lr()) << "]" << endl;
		os << "" << endl;
		//Gamma
		os << "Statistics of Call Option: Monte Carlo - Gamma(LR-LR):" << endl;
		os << "Estimate of Call Option Gamma: " << gamma_mc_lr_lr << endl;
		os << "Actual Call Option Gamma " << gamma_bs() << endl;
		os << "Bias of estimator: " << gamma_bs() - gamma_mc_lr_lr << endl;
		os << "Variance of estimator: " << variance_gamma_lr_lr() << endl;
		os << "Mean Square Error: " << variance_gamma_lr_lr() + pow((gamma_bs() - gamma_mc_lr_lr), 2) << endl;
		os << "99.7% Confidence Interval [" << gamma_mc_lr_lr - 3 * sqrt(variance_gamma_lr_lr()) << ", " << gamma_mc_lr_lr + 3 * sqrt(variance_gamma_lr_lr()) << "]" << endl;
		os << "" << endl;

		os << "Statistics of Call Option: Monte Carlo - Gamma(LR-PW):" << endl;
		os << "Estimate of Call Option Gamma: " << gamma_mc_lr_pw << endl;
		os << "Actual Call Option Gamma " << gamma_bs() << endl;
		os << "Bias of estimator: " << gamma_bs() - gamma_mc_lr_pw << endl;
		os << "Variance of estimator: " << variance_gamma_lr_pw() << endl;
		os << "Mean Square Error: " << variance_gamma_lr_pw() + pow((gamma_bs() - gamma_mc_lr_pw), 2) << endl;
		os << "99.7% Confidence Interval [" << gamma_mc_lr_pw - 3 * sqrt(variance_gamma_lr_pw()) << ", " << gamma_mc_lr_pw + 3 * sqrt(variance_gamma_lr_pw()) << "]" << endl;
		os << "" << endl;

		os << "Statistics of Call Option: Monte Carlo - Gamma(PW-LR):" << endl;
		os << "Estimate of Call Option Gamma: " << gamma_mc_pw_lr << endl;
		os << "Actual Call Option Gamma " << gamma_bs() << endl;
		os << "Bias of estimator: " << gamma_bs() - gamma_mc_pw_lr << endl;
		os << "Variance of estimator: " << variance_gamma_pw_lr() << endl;
		os << "Mean Square Error: " << variance_gamma_pw_lr() + pow((gamma_bs() - gamma_mc_pw_lr), 2) << endl;
		os << "99.7% Confidence Interval [" << gamma_mc_pw_lr - 3 * sqrt(variance_gamma_pw_lr()) << ", " << gamma_mc_pw_lr + 3 * sqrt(variance_gamma_pw_lr()) << "]" << endl;
		os << "" << endl;

		return os;
	}
};


////////////////////////////////////////////////////////////////////////////////////////////////
//				Distribution of S_n : Heterogeneous Portfolio								  //
////////////////////////////////////////////////////////////////////////////////////////////////

// With functions
// With C++ tool : to generate standard Uniform random variable
//default_random_engine generator;
//uniform_real_distribution<double> distribution(0.0, 1.0);


///////////////////////////////// Andersen Sidenius Basu algo


double Andersen_sidenius_Basu(int k, unsigned int j, int n) {
	/* this algorithm takes 3 integers as arguments (k,j and n) and compute the distribution P(S_n = k) through the recursive relation of the Andersen Sidenius Basu
	Algorithm given by : P(S_{j}=k)P(L_{j+1}=0) + P(S_{j}=k-1)P(L_{j+1}=1). This algorithm is a recursive one.*/

	//Initial condition
	if ((k > n) | (j < 0) | (k < 0)) {
		return 0;
	}
	else {
		if ((j == 0) & (k == 0)) {
			return 1;
		}
		else {
			//in case k=j, then it means we have k successes so we can compute faster the distribution
			if (k == j) {
				double Pk = 1;
				for (unsigned int i = 0; i <= j; i++) {
					double P_L = 1 - ((double)i / ((double)n + 1));
					Pk *= P_L;
				}
				return Pk;
			}
			else {
				//Recursive relation
				double P_L = 1 - ((double)j / ((double)n + 1));
				return Andersen_sidenius_Basu(k, j - 1, n) * (1 - P_L) + Andersen_sidenius_Basu(k - 1, j - 1, n) * (P_L);
			}
		}
	}
}




///////////////////////////////// Hull White algorithm

double Hull_White(unsigned int n, unsigned int k) {
	/*this algorithm takes 2 integers as arguments (n and k) and compute the distribution P(S_n = k) through the Hull and White Algorithm*/
	//Computation of pi
	double pi_n = 1;
	for (unsigned int i = 1; i <= n; i++) {
		double PL = 1 - ((double)i / ((double)n + 1));
		pi_n *= (1 - PL);
	}
	// Case k=0
	if (k == 0) {
		return pi_n;
	}
	// Case k>0
	double cj;
	vector<double> Vn;
	vector<double> Un;
	// Recursion
	for (unsigned int i = 0; i < k; i++) {
		//TEST cout << "k = " << i << endl;

		// Computation of Vi
		double Vi = 0;
		for (unsigned int j = 1; j <= n; j++) {
			double PL = 1 - ((double)j / ((double)n + 1));
			cj = PL / (1 - PL);
			//TEST cout << "P = " << PL << endl;
			//TEST cout << "ci = " << ci << endl;
			Vi += pow(cj, i + 1);
		}
		//TEST cout << "i = " << i << "; Vi =" << Vi << endl;
		Vn.push_back(Vi);

		// Computation of Ui
		double Ui = 0;
		if (Un.size() == 0) {
			Un.push_back(1);
		}

		for (unsigned int l = 0; l < Vn.size(); l++) {
			//TEST cout << "V" << l << " = " << Vn[l]<<endl;
			//TEST cout << "U" << Un.size() - l-1 << " = " << Un[Un.size() - l-1] << endl;
			Ui += pow(-1, l) * (Vn[l] * Un[Un.size() - l - 1]);
		}
		Ui = Ui / (Un.size());
		Un.push_back(Ui);
	}

	return pi_n * Un[Un.size() - 1];
}

///////////////////////////////// Monte Carlo Simulation


class het_probabilities {
	/*This class allows to simulate the sum of random variables Sn according to the probabilities pi and compute its distribution.
	It also presents the functions to compute the distribution of Sn via Andersen Sidenius Basu algorithm and the relative error*/
	int simulations;
	unsigned int max_k;
	vector<double> probs;
	vector<double> and_sid_bas_probs;
	vector<double> hull_white_probs;
	vector<double> error_asd_dist;
	vector<double> error_hw_dist;
public:
	het_probabilities(const int& sims = 1, const int& k = 1) {
		simulations = sims;
		max_k = k;
		vector<double> temp(max_k + 1, 0);
		probs = temp;
		and_sid_bas_probs = temp;
		hull_white_probs = temp;
		error_hw_dist = temp;
		error_asd_dist = temp;
	}

	void simulation_probs() {
		start_time = chrono::system_clock::now();
		random_gen U(simulations * max_k);
		for (int sims = 1; sims <= simulations; sims++) {
			int S_n = 0;
			for (unsigned int i = 1; i <= max_k; i++) {
				double rand_n = U.get_unif()[i*sims - 1];
				double p_i = 1 - (double(i) / (max_k + 1));
				if (rand_n < p_i) {
					S_n += 1;
				}
			}
			probs[S_n] += 1;
		}
		for (unsigned int j = 0; j <= max_k; j++) {
			probs[j] = double(probs[j]) / simulations;
		}
		double time_taken = run_time_sec();
		cout << "Run time for Monte Carlo Simulation for " << simulations << " simulations: " << time_taken << " seconds" << endl;
	}
	void And_Sid_Bas() {
		start_time = chrono::system_clock::now();
		for (unsigned int k = 0; k <= max_k; k++) {
			double p = Andersen_sidenius_Basu(k, max_k, max_k);
			and_sid_bas_probs[k] = (p);
		}
		double time_taken = run_time_sec();
		cout << "Run time for Andersen Sidenius Basu Algorithm: " << time_taken << " seconds" << endl;
	}
	void Hull_White_dis() {
		start_time = chrono::system_clock::now();
		for (unsigned int k = 0; k <= max_k; k++) {
			double p = Hull_White(max_k, k);
			hull_white_probs[k] = (p);
		}
		double time_taken = run_time_sec();
		cout << "Run time for Hull White Algorithm: " << time_taken << " seconds" << endl;
	}
	void errors_ASD_computation() {
		for (unsigned int k = 0; k <= max_k; k++) {
			double e = abs(probs[k] - and_sid_bas_probs[k]) / and_sid_bas_probs[k];
			error_asd_dist[k] = e;
		}
	}
	void errors_HW_computation() {
		for (unsigned int k = 0; k <= max_k; k++) {
			double e = abs(probs[k] - hull_white_probs[k]) / hull_white_probs[k];
			error_hw_dist[k] = e;
		}
	}
	void print() {
		cout << "\n\n-------------------------------------------PROBABILITIES OF S_n----------------------------------------------" << endl;
		cout << "MONTE CARLO (MC) PROBABILITIES: \n";
		for (unsigned int i = 0; i <= max_k; i++) {
			cout << "Simulation result : P(S" << max_k << "=" << i << ") = " << probs[i] << endl;
		}
		cout << "\n";
		cout << "-------------------------------------------------------------------------------------------------------------" << endl;
		cout << "ANDERSEN SIDENIUS BASU PROBABILITIES: \n";
		for (unsigned int i = 0; i <= max_k; i++) {
			cout << "Andersen Sidenius Basu Algorithm result : P(S" << max_k << "=" << i << ") = " << and_sid_bas_probs[i] << endl;
		}
		cout << "\n";
		cout << "-------------------------------------------------------------------------------------------------------------" << endl;
		cout << "HULL WHITE PROBABILITIES: \n";
		for (unsigned int i = 0; i <= max_k; i++) {
			cout << "Hull White Algorithm result : P(S" << max_k << "=" << i << ") = " << hull_white_probs[i] << endl;
		}
		cout << "\n";
		cout << "-------------------------------------------------------------------------------------------------------------" << endl;
		cout << "MC vs. ANDERSEN SIDENIUS BASU RELATIVE ERRORS: \n";
		for (unsigned int i = 0; i <= max_k; i++) {
			cout << "Andersen Sidenius Basu Algorithm relative error for : P(S" << max_k << "=" << i << ") is " << error_asd_dist[i] << endl;
		}
		cout << "\n";
		cout << "-------------------------------------------------------------------------------------------------------------" << endl;
		cout << "MC vs. HULL WHITE RELATIVE ERRORS: \n";
		for (unsigned int i = 0; i <= max_k; i++) {
			cout << "Hull White Algorithm relative error for : P(S" << max_k << "=" << i << ") is " << error_hw_dist[i] << endl;
		}
		cout << "\n";
		cout << "-------------------------------------------------------------------------------------------------------------" << endl;
	}

	void write(string filename) {
		/*string filename;
		cin >> filename;*/
		ofstream myfile;
		myfile.open(filename + ".csv");
		myfile << "k;Simulation;Andersen Sidenius Basu Algorithm result;Andersen Sidenius Basu Algorithm error;Hull White Algorithm result;Hull White Algorithm error .\n";
		for (unsigned int i = 0; i <= max_k; i++) {
			myfile << i << ";" << probs[i] << ";" << and_sid_bas_probs[i] << ";" << error_asd_dist[i] << ";" << hull_white_probs[i] << ";" << error_hw_dist[i] << ".\n";
		}
	}
};



////////////////////////////////////////////////////////////////////////////////////////////////
//								Tails Distribution study									  //
////////////////////////////////////////////////////////////////////////////////////////////////



class tail {
protected:
	unsigned int n, N;
public:
	tail(int nb_bond, int nb_sims) {
		n = nb_bond;
		N = nb_sims;
	}
	virtual int generate_Sn() = 0;
	virtual double distribution(double x) = 0;
	virtual double Cramer(double x) = 0;
	virtual double CLT(double x) = 0;
	virtual double gamma(double x) = 0;
	void print(double xc, double xclt) {
		cout << "n = " << n << " ; " << endl;
		cout << "Distribution : P(Sn > n*x) = " << distribution(xclt) << endl;
		cout << "--------- Verification of Cramer Theorem with x = " << xc << "---------" << endl;
		cout << "Cramer Theorem : ln(P(Sn>n*x))/n = " << Cramer(xc) << endl;
		cout << "-Gamma*(x) = " << -gamma(xc) << endl;
		cout << "--------- Verification of Central Limit Theorem with x = " << xclt << " ---------" << endl;
		cout << "ln(P((Sn - n*mu)/(sqrt(n)sigma)>x))/n = " << CLT(xclt) << endl;
		cout << "1-NormalCDF(x) =" << 1 - normalCDF(xclt) << endl;
	}
	void write(double xc, double xclt, ofstream& myfile0) {
		myfile0 << n << ";" << xc << ";" << distribution(xc) << ";" << Cramer(xc) << ";" << -gamma(xc) << ";" << xclt << ";" << CLT(xclt) << ";" << 1 - normalCDF(xclt) << ".\n";
	}
};
///////////////////////////////// Bernoulli Distribution ///////////////////////////////////////// 

///////////////////////////////// Monte Carlo simulation 
// Poisson generator by C++ tool
//default_random_engine generator_bernoulli;
// Poisson generator by C++ tool
//default_random_engine generator_poisson;


class bernoulli : public tail {
	double p;
public:
	bernoulli(int nb_bond, int nb_sims, double parameter) :tail(nb_bond, nb_sims) {
		p = parameter;
	}
	int generate_Sn() {
		/* This function takes a double p and an integer n as arguments and return the sum of n Bernoulli random variables of parameter p*/
	// Tool to generate the Bernoulli random  variable of parameter p
		bernoulli_distribution bernoulli_distribution(p);
		// Initialisation of the sum
		int S_n = 0;
		// Generation the random variables and computation of the sum
		//// Initial code
		/*for (unsigned int i = 0; i < n; i++) {
		//	int Xi_Ber = bernoulli_distribution(generator_bernoulli);
		//	//TEST cout << Xi_poi << endl;
		//	S_n += Xi_Ber;
		//}*/
		random_gen U(n);
		vector<int> Bernoulli = U.bernoulli(p);
		for (unsigned int i = 0; i < Bernoulli.size(); i++) {
			
			S_n += Bernoulli[i];
		}
		//cout << "Sn = " << S_n << endl;
		return S_n;
	}
	double distribution(double x) {
		/*This function takes 2 doubles x and p and 2 integers N and n as arguments and returns the estimated value of P(S_n > nx) with N simulations*/

	// Initialisation : we store the values of the sum S_n in a vector
		vector<int> Sn_values;
		// Simulation of S_n, N times
		for (unsigned int i = 0; i < N; i++) {
			//cout << "x=" << x << endl;
			int S_n = generate_Sn();
			Sn_values.push_back(S_n);
		}
		// Count how many times S_n has been superior to nx
		int count = 0;
		for (unsigned int k = 0; k < Sn_values.size(); k++) {
			if (Sn_values[k] >= n * x) {
				count += 1;
			}
		}
		return (double)count / (double)N;
	}
	double Cramer(double x) {
		/*This function takes 2 doubles x and p and 2 integers N and n as arguments and returns ln(P(S_n > nx))/n*/
		return log(distribution(x)) / (double)n;
	}
	double gamma(double x) {
		/*this function takes 2 doubles as arguments x and lambda and returns the value Gamma*(x) for a Poisson distribution of parameter lambda*/
		double inf = numeric_limits<double>::infinity();
		if (x < 0) {
			return inf;
		}
		else {
			if (x == 0) {
				return -log(1 - p);
			}
			else {
				if ((x > 0)& (x < 1)) {
					return x * log(x * (1 - p) / (p * (1 - x))) - log(x * (1 - p) / (1 - x) + 1 - p);
				}
				else {
					if (x == 1) {
						return -log(p);
					}
					else {
						if (x > 1) {
							return inf;
						}
					}
				}
			}
		}
	}
	double CLT(double x) {
		/*This function takes 3 doubles x, N and lambda and 1 integer n as arguments and returns the following probability P(>)*/
		double X = (p * (1 - p) * x) / sqrt((double)n) + p;
		//cout << X << endl;
		return distribution(X);
	}
};


///////////////////////////////// Poisson Distribution ///////////////////////////////////////// 

class poisson : public tail {
	double lambda;
public:
	poisson(int nb_bond, int nb_sims, double parameter) :tail(nb_bond, nb_sims) {
		lambda = parameter;
	}
	int generate_Sn() {
		/* This function takes a double lambda and an integer n as arguments and return the sum of n poisson random variables of parameter lambda*/
	// Tool to generate the poisson random  variable of parameter lambda
		//poisson_distribution<int> poisson_distribution(lambda);
		// Initialisation of the sum
		int S_n = 0;
		// Generation the random variables and computation of the sum
		// Initial code
		/*for (unsigned int i = 0; i < n; i++) {
			int Xi_poi = poisson_distribution(generator_poisson);
			//TEST cout << Xi_poi << endl;
			S_n += Xi_poi;
		}*/
		random_gen U(n);
		vector<int> Poisson = U.poisson(lambda);
		for (unsigned int i = 0; i < Poisson.size(); i++) {
			S_n += Poisson[i];
		}
		return S_n;
	}
	double distribution(double x) {
		/*This function takes 2 doubles x and lambda and 2 integers N and n as arguments and returns the estimated value of P(S_n > nx) with N simulations*/

	// Initialisation : we store the values of the sum S_n in a vector
		vector<int> Sn_values;
		// Simulation of S_n, N times
		for (unsigned int i = 0; i < N; i++) {
			int S_n = generate_Sn();
			Sn_values.push_back(S_n);
		}
		// Count how many times S_n has been superior to nx
		int count = 0;
		for (unsigned int k = 0; k < Sn_values.size(); k++) {
			if (Sn_values[k] >= n * x) {
				count += 1;
			}
		}
		return (double)count / (double)N;
	}
	double Cramer(double x) {
		/*This function takes 2 doubles x and lambda and 2 integers N and n as arguments and returns ln(P(S_n > nx))/n*/
		return log(distribution(x)) / (double)n;
	}
	double gamma(double x) {
		/*this function takes 2 doubles as arguments x and lambda and returns the value Gamma*(x) for a Poisson distribution of parameter lambda*/
		return log(x / lambda) * x - x + lambda;
	}
	double CLT(double x) {
		/*This function takes 3 doubles x, N and lambda and 1 integer n as arguments and returns the following probability P(>)*/
		double X = (sqrt(lambda) * x) / sqrt((double)n) + lambda;
		//cout << X << endl;
		return distribution(X);
	}
};







////////////////////////////////////////////////////////////////////////////////////////////////
//									   Main Functions    									  //
////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////// Main Function for the global task

int main0() {
	double S0, K, T, R, sigma;
	int choice_m0 = -1;
	while (choice_m0 != 99) {
		cout << "\n\n--------------------------------------------Global task----------------------------------------------" << endl;
		cout << "Welcome to the Global Task of the Coursework section. \nHere we compute the price of a call option and its greeks using:";
		cout << "\n----------- (1) Monte-Carlo, (2) Black - Scholes formula ----------\n" << endl;
		cout << "\nSeveral statistics are provided like the mean, the variance, the execution time...\n" << endl;
		cout << "(*)  Continue (any number, excl. 99); \n(99) Return to previous menu\n";
		cout << "Choice: ";
		cin >> choice_m0;
		if (choice_m0 != 99) {
			cout << "\nPlease choose a value for the initial price S0 of the call option" << endl;
			cout << "Choice for 'S' (100 for example): ";
			cin >> S0;
			cout << "\nPlease choose a value for the strike K of the call option" << endl;
			cout << "Choice for 'K' (110 for example):";
			cin >> K;
			cout << "\nPlease choose a value for the interest rate rho of the call option" << endl;
			cout << "Choice for 'rho'(0.05 for example): ";
			cin >> R;
			cout << "\nPlease choose a value for the maturity T of the call option" << endl;
			cout << "Choice for 'T' in number of ??? (??? for example):";
			cin >> T;
			cout << "\nPlease choose a value for the volatility sigma of the call option" << endl;
			cout << "Choice for 'sigma' (0.2 for example): ";
			cin >> sigma;
			cout << "\nStarting simulation/ calculations..." << endl;

			ofstream ofile("price.txt"); // creates an ofstream called ofile
			if (!ofile) {
				cout << "error opening file";
				exit(1); // error opening file
			}
			Call_EU s(S0, K, R, T, sigma);
			s.MonteCarlo(1000000);
			s.Statistics(cout);
			s.Statistics(ofile);
		}
		else if (choice_m0 == 99) {
			return 0;
		}
		else {
			cout << "\n\n" << endl;
			cout << "///////////////////////////////////////" << endl;
			cout << "// Invalid choice, please try again. //" << endl;
			cout << "///////////////////////////////////////\n\n\n\n";
		}



	}
}


//////////////////////////////////// Main Function for the distribution of Sn

int main1()
{
	//////////////
	int sims, blabla;

	int choice_m1 = -1; // safeguard to prevent any unwanted behaviour in program
	while (choice_m1 != 99) {
		cout << "\n\n--------------------------------------------DISTRIBUTION OF S_n----------------------------------------------" << endl;
		cout << "Welcome to the sum of random variables (S_n) distribution section. \nHere we calculate the probability distributions for the sum of random variables using:";
		cout << "\n----------- (1) Monte-Carlo, (2) ANDERSEN SIDENIUS BASU Algorithm and the (3) HULL WHITE Algorithm ----------\n" << endl;
		cout << "(*)  Continue (any number, excl. 99); \n(99) Return to previous menu\n";
		cout << "Choice: ";
		cin >> choice_m1;

		if (choice_m1 != 99) {
			cout << "\nPlease choose an integer value 'n' for the value S_n = L1+L2+...Ln" << endl;
			cout << "Choice (for 'n' > 0): ";
			cin >> blabla;
			cout << "\nNumber of simulations: ";
			cin >> sims;

			cout << "\nStarting simulation/ calculations..." << endl;

			// Compute Sn
			het_probabilities probabilities(sims, blabla);
			// Distribution
			probabilities.simulation_probs();
			// Andersen Sidenius Basu Algorithm
			probabilities.And_Sid_Bas();
			// Hull White Algorithm
			probabilities.Hull_White_dis();

			//Errors computation
			probabilities.errors_ASD_computation();
			probabilities.errors_HW_computation();

			//Print
			probabilities.print();
			// Write
			//auto name = to_string(sims);
			//probabilities.write(name + "simulations");
			cout << "\nSimulation/ calculations complete." << endl;
		}
		else if (choice_m1 == 99) {
			return 0;
		}
		else {
			cout << "\n\n" << endl;
			cout << "///////////////////////////////////////" << endl;
			cout << "// Invalid choice, please try again. //" << endl;
			cout << "///////////////////////////////////////\n\n\n\n";
		}

	}

	return 0;
}


//////////////////////////////////// Main Function for the tail of Sn

int main2() {
	double x = -2;
	vector<double> x_clt;
	for (unsigned int i = 0; i <= 40; i++) {
		x += 0.1;
		x_clt.push_back(x);
	}

	int choice_m2 = -1; // safeguard to prevent any unwanted behaviour in program
	while (choice_m2 != 99) {
		int n, N = 1000;// n will go to infinity and N is the number of simulations
		cout << "\n\n---------------------------------------------TAIL DISTRIBUTIONS----------------------------------------------" << endl;
		cout << "Welcome to the tail distribution section, where we study properties of tail distributions. \nTwo examples can be studied: Please choose the one you want: \n(0)  Bernoulli distribution; \n(1)  Poisson distribution; \n(99) Return to previous menu;" << endl;
		cout << "Choice: ";
		cin >> choice_m2;
		if (choice_m2 == 0) {
			// Bernoulli
			// tails parameters : 
			vector<double> x_bernoulli_cramer;
			cout << "\n" << endl;
			cout << "////////////////////////////////" << endl;
			cout << "//// BERNOULLI distribution ////" << endl;
			cout << "////////////////////////////////" << endl;
			cout << "\n" << endl;
			cout << "Please choose the probability p of a Bernoulli distribution." << endl;
			cout << "10 iterations will be run using 10,000 values greater than or equal to E[x]" << endl;
			cout << "(Results will be output to a csv. NOTE: Previous csv will be overwritten)" << endl;
			cout << "Probability (between 0 and 1): ";
			double p;
			cin >> p;
			cout << "\nStarting simulations..." << endl;
			double xb = p;
			ofstream myfile0;
			myfile0.open("Bernoulli_Tail_distribution.csv");
			myfile0 << "n;x;Pber(Sn>x);Cramer;-Gamma*;x;Central_limit;1-NormalCDF.\n";
			for (unsigned int i = 0; i <= 40; i++) {
				xb += (1-p)/41;
				//cout << xb << endl;
				x_bernoulli_cramer.push_back(xb);
			}
			for (int i = 1; i < 10; i++) {
				n = i * 100;
				bernoulli B(n, N, p);
				for (unsigned int j = 0; j < x_bernoulli_cramer.size(); j++) {
					double x_bc = x_bernoulli_cramer[j];
					double x_clt_double = x_clt[j];
					//cout << x_bc << endl;
					//bernoulli B(n, N, p);
					/*if (j % 100 == 0) {
						cout << "\n";
						cout << "-------------------------------------------------------------------------------------------------------------" << endl;
						cout << "Iteration:" << i << endl;
						cout << "-------------------------------------------------------------------------------------------------------------" << endl;
						B.print(x_bc, x_clt_double);
					}*/
					cout << "\n";
					cout << "-------------------------------------------------------------------------------------------------------------" << endl;
					cout << "Iteration:" << i << endl;
					cout << "-------------------------------------------------------------------------------------------------------------" << endl;
					B.print(x_bc, x_clt_double);
					B.write(x_bc, x_clt_double, myfile0);
				}

				//P.write(x, myfile1);
			}
			myfile0.close();
			cout << "\nSimulation complete. \nPlease find output in the file: Bernoulli_Tail_distribution.csv\n" << endl;

		}
		else if (choice_m2 == 1) {
			// Poisson
			// Parameters : 
			N = 10000;
			vector<double> x_poisson_cramer;
			cout << "\n" << endl;
			cout << "////////////////////////////////" << endl;
			cout << "///// POISSON distribution /////" << endl;
			cout << "////////////////////////////////" << endl;
			cout << "\n" << endl;
			cout << "Please choose the parameter lambda for the Poisson distribution.\n";
			cout << "10 iterations will be run using 10,000 values greater than or equal to E[x]" << endl;
			cout << "(Results will be output to a csv. NOTE: Previous csv will be overwritten)" << endl;
			cout << "Lambda: ";
			double lambda;
			cin >> lambda;
			cout << "\nStarting simulations..." << endl;
			double xp = lambda;
			ofstream myfile1;
			myfile1.open("Poisson_Tail_distribution.csv");
			myfile1 << "n;Pber(Sn>x);Cramer;Gamma*;Central_limit;1-Normal.\n";
			for (unsigned int i = 0; i <= 40; i++) {
				xp += 0.1;
				x_poisson_cramer.push_back(xp);
			}
			for (int i = 1; i < 10; i++) {
				n = i * 10;
				poisson P(n, N, lambda);
				for (unsigned int j = 0; j < x_poisson_cramer.size(); j++) {
					double x_pc = x_poisson_cramer[j];
					double x_clt_double = x_clt[j];
					//poisson P(n, N, lambda);
					/*if (j % 10 == 0) {
						cout << "\n";
						cout << "-------------------------------------------------------------------------------------------------------------" << endl;
						cout << "Iteration:" << i << endl;
						cout << "-------------------------------------------------------------------------------------------------------------" << endl;
						P.print(x_pc, x_clt_double);
					}*/
					cout << "\n";
					cout << "-------------------------------------------------------------------------------------------------------------" << endl;
					cout << "Iteration:" << i << endl;
					cout << "-------------------------------------------------------------------------------------------------------------" << endl;
					P.print(x_pc, x_clt_double);
					P.write(x_pc, x_clt_double, myfile1);
				}


			}
			myfile1.close();
			cout << "\nSimulation complete. \nPlease find output in the file: Poisson_Tail_distribution.csv\n" << endl;
		}
		else if (choice_m2 == 99) {
			return 0;
		}
		else {
			cout << "\n\n" << endl;
			cout << "///////////////////////////////////////" << endl;
			cout << "// Invalid choice, please try again. //" << endl;
			cout << "///////////////////////////////////////\n\n\n\n";
		}
	}

	return 0;

}

//int main2() {
//	int n, N = 10000;// n will go to infinity and N is the number of simulations
////			// Bernoulli
////			// tails parameters : 
//	double x = -2;
//	vector<double> x_clt;
//	for (unsigned int i = 0; i <= 40; i++) {
//				x += 0.1;
//				x_clt.push_back(x);
//			}
//	vector<double> x_bernoulli_cramer;
//	double p=0.3;
//	double xb = p;
//	for (unsigned int i = 0; i <= 40; i++) {
//		xb += (1-p)/41;
//		//cout << xb << endl;
//		x_bernoulli_cramer.push_back(xb);
//		}
//	for (int i = 1; i < 10; i++) {
//		n = i * 100;
//		bernoulli B(n, N, p);
//		for (unsigned int j = 0; j < x_bernoulli_cramer.size(); j++) {
//			double x_bc = x_bernoulli_cramer[j];
//			double x_clt_double = x_clt[j];
//			cout << x_bc << endl;
//			B.print(x_bc, x_clt_double);
//		}
//	}
//	return 0;
//}

//////////////////////////////////// Main Function

int main() {
	int choice_m = -1;
	while (choice_m != 99) {
		cout << "\n\n--------------------------------------------------MAIN MENU--------------------------------------------------" << endl;
		cout << "Please enter a choice for the SIMULATION METHOD: \n(0)  Global task; \n(1)  Distribution of Sn; \n(2)  Tails of distribution; \n(99) Exit program;" << endl;
		cout << "Choice: ";
		cin >> choice_m;
		if (choice_m == 0) {
			main0();
			cout << "UNDER CONSTRUCTION... SAMI/ HITCHAM PLEASE SEND :P\n";
		}
		else if (choice_m == 1) {
			main1();
		}
		else if (choice_m == 2) {
			main2();
		}
		else if (choice_m == 99) {
			cout << "\n\n" << endl;
			cout << "////////////////////////////////" << endl;
			cout << "////// EXITING. Good bye. //////" << endl;
			cout << "////////////////////////////////" << endl;
			cout << "\n\n" << endl;
		}
		else {
			cout << "\n\n" << endl;
			cout << "///////////////////////////////////////" << endl;
			cout << "// Invalid choice, please try again. //" << endl;
			cout << "///////////////////////////////////////\n\n\n\n";
		}
	}

	return 0;
}

