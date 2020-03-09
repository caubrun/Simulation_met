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
#include <iomanip>
#include <functional> 
#include <fstream>
#include <string>

#include <ctime>
#include <ratio>
#include <chrono>

using namespace std;

////////////////////////////////////////////////////////////////////////////////////////////////
//									  Timer functions 										  //
////////////////////////////////////////////////////////////////////////////////////////////////


// GLOBAL VARIABLE TO STORE STARTING TIME:
chrono::high_resolution_clock::time_point start_time;

void start_timing() {
	start_time = chrono::high_resolution_clock::now();
}

double run_time_sec() {
	chrono::high_resolution_clock::time_point t2 = chrono::high_resolution_clock::now();
	chrono::duration<double> time_span = chrono::duration_cast<std::chrono::microseconds>(t2 - start_time);
		//chrono::duration_cast<chrono::duration<double>>(t2 - start_time);
	double run_duration = time_span.count()/1000000;
	return run_duration;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//				Distribution of S_n : Heterogeneous Portfolio								  //
////////////////////////////////////////////////////////////////////////////////////////////////

// With functions
// With C++ tool : to generate standard Uniform random variable
default_random_engine generator;
uniform_real_distribution<double> distribution(0.0, 1.0);

///////////////////////////////// Andersen Sidenius Basu algo


/* This looks like it clogs up a lot of memory, as it creates multiple instances of the functions, goes through the same checks again, etc.... rather than  completing the function in a loop*/
/* I cann suggest an improvement, but maybe in the end... let's focus on other stuff for now*/
/* NOTE: I THINK THERE IS AN ERROR BELOW. */
/* Ok, I misunderstood what that part of the code was doing*/

double Andersen_sidenius_Basu(int k, int j, int n) {
	/* this algorithm takes 3 integers as arguments (k,j and n) and compute the distribution P(S_n = k) through the recursive relation of the Andersen Sidenius Basu
	Algorithm given by : P(S_{j}=k)P(L_{j+1}=0) + P(S_{j}=k-1)P(L_{j+1}=1). This algorithm is a recursive one.*/

	//Initial condition
	if ((k > n) | (j < 0) | (k < 0)) {
		//TEST cout << "0 case : n=" << n << " and k = " << k << endl;
		return 0;
	}
	else {
		if ((j == 0) & (k == 0)) {
			//TEST cout << "1 case" << endl;
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
				//TEST cout << P_L << "k = " << k << " ; n =" << n << endl;
				return Andersen_sidenius_Basu(k, j - 1, n) * (1 - P_L) + Andersen_sidenius_Basu(k - 1, j - 1, n) * (P_L);
			}
		}
	}
}




///////////////////////////////// Hull White algorithm

double Hull_White(int n, int k) {
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

	/*for (unsigned int i = 1; i <= n; i++) {
		double PL = 1 - ((double)i / ((double)n + 1));
		pi_n *= (1 - PL);
	}*/
	return pi_n * Un[Un.size() - 1];
}

///////////////////////////////// Monte Carlo Simulation

// With a class
class het_probabilities {
	int simulations;
	int max_k;
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

	//void simulation_probs(vector<double>& sim_probs) {
	void simulation_probs() {
		//vector<double> sim_outcome;
		for (int sims = 1; sims <= simulations; sims++) {
			int S_n = 0;
			for (int i = 1; i <= max_k + 1; i++) {
				random_device rd;  //Will be used to obtain a seed for the random number engine
				mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
				uniform_real_distribution<> dis(0.0, 1.0);
				// Use dis to transform the random unsigned int generated by gen into a 
				// double in [1, 2). Each call to dis(gen) generates a new random double
				double rand_n = dis(gen);
				double p_i = 1 - (double(i) / (max_k + 1));
				//cout << rand_n << " " << p_i << endl;
				if (rand_n < p_i) {
					S_n += 1;
				}
			}
			//sim_outcome[sims] = S_n;
			//cout << "S_n=" << S_n << endl;
			probs[S_n] += 1;
			//cout << "probs[S_n]=" << probs[S_n] << endl;
		}
		//transform(sim_probs.begin(), sim_probs.end(), sim_probs.begin(), bind(divides<int>(simulations), placeholders::_1, 3));
		for (int j = 0; j <= max_k; j++) {
			probs[j] = double(probs[j]) / simulations;
		}
	}
	void And_Sid_Bas() {
		start_time;
		for (unsigned int k = 0; k <= max_k; k++) {
			double p = Andersen_sidenius_Basu(k, max_k, max_k);
			and_sid_bas_probs[k] = (p);
		}
		double time_taken = run_time_sec();
		cout << "Run time:" << time_taken << " seconds" << endl;

	}
	void Hull_White_dis() {
		for (unsigned int k = 0; k <= max_k; k++) {
			double p = Hull_White(max_k, k);
			hull_white_probs[k] = (p);
		}
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
		for (int i = 0; i <= max_k; i++) {
			cout << "Simulation result : P(S" << max_k << "=" << i << ") = " << probs[i] << endl;
			cout << "Andersen Sidenius Basu Algorithm result : P(S" << max_k << "=" << i << ") = " << and_sid_bas_probs[i] << endl;
			cout << "Andersen Sidenius Basu Algorithm relative error for : P(S" << max_k << "=" << i << ") is " << error_asd_dist[i] << endl;
			cout << "Hull White Algorithm result : P(S" << max_k << "=" << i << ") = " << hull_white_probs[i] << endl;
			cout << "Hull White Algorithm relative error for : P(S" << max_k << "=" << i << ") is " << error_hw_dist[i] << endl;
		}
	}

	void write(string filename) {
		/*string filename;
		cin >> filename;*/
		ofstream myfile;
		myfile.open(filename + ".csv");
		myfile << "k;Simulation;Andersen Sidenius Basu Algorithm result;Andersen Sidenius Basu Algorithm error;Hull White Algorithm result;Hull White Algorithm error .\n";
		for (int i = 0; i <= max_k; i++) {
			myfile << i << ";" << probs[i] << ";" << and_sid_bas_probs[i] << ";" << error_asd_dist[i] << ";" << hull_white_probs[i] << ";" << error_hw_dist[i] << ".\n";
		}
	}
};



////////////////////////////////////////////////////////////////////////////////////////////////
//								Tails Distribution study									  //
////////////////////////////////////////////////////////////////////////////////////////////////
double normalCDF(double x)
{/*normal cdf*/
	return erfc(-x / sqrt(2)) / 2;
}


class tail {
protected:
	int n, N;
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
	void write(double x, ofstream &myfile0) {
		myfile0 << "n;Pber(Sn>x);Cramer;Gamma*;Central_limit;1-Normal.\n";
		myfile0 << n << ";" << distribution(x) << ";" << Cramer(x) << ";" << -gamma(x) << ";" << CLT(x) << ";" << normalCDF(x) << ".\n";
	}
};
///////////////////////////////// Bernoulli Distribution ///////////////////////////////////////// 

///////////////////////////////// Monte Carlo simulation 
// Poisson generator by C++ tool
default_random_engine generator_bernoulli;
// Poisson generator by C++ tool
default_random_engine generator_poisson;


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
		for (unsigned int i = 0; i < n; i++) {
			int Xi_Ber = bernoulli_distribution(generator_bernoulli);
			//TEST cout << Xi_poi << endl;
			S_n += Xi_Ber;
		}
		return S_n;
	}
	double distribution(double x) {
		/*This function takes 2 doubles x and p and 2 integers N and n as arguments and returns the estimated value of P(S_n > nx) with N simulations*/

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
		poisson_distribution<int> poisson_distribution(lambda);
		// Initialisation of the sum
		int S_n = 0;
		// Generation the random variables and computation of the sum
		for (unsigned int i = 0; i < n; i++) {
			int Xi_poi = poisson_distribution(generator_poisson);
			//TEST cout << Xi_poi << endl;
			S_n += Xi_poi;
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
//									   Main Function    									  //
////////////////////////////////////////////////////////////////////////////////////////////////


int main()
{
	//////////////
	int sims, blabla = 10;

	for (unsigned int i = 1; i <= 10; i++) {
		sims = i * 1000;
		het_probabilities probabilities(sims, blabla);
		probabilities.simulation_probs();
		probabilities.And_Sid_Bas();
		probabilities.Hull_White_dis();
		probabilities.errors_ASD_computation();
		probabilities.errors_HW_computation();
		//probabilities.print();
		auto name = to_string(sims);
		probabilities.write(name + "simulations");
	}


	////////////////
	// Tails of Distribution of Sn

	int N = 100; // number of Simulation
	int n = 10; // n

	// tails proprieties
	double x = 5.02;
	double lambda = 5;
	double xp = 0.51;
	double p = 0.5;

	ofstream myfile0;
	myfile0.open("Bernoulli_Tail_distribution.csv");
	ofstream myfile1;
	myfile1.open("Poisson_Tail_distribution.csv");
	for (int i = 1; i < 1000; i++) {
		n = i * 1000;
		bernoulli B(n, N, p);
		B.write(xp, myfile0);
		poisson P(n, N, lambda);
		P.write(x, myfile1);
	}
	myfile0.close();
	myfile1.close();

	return 0;

}

