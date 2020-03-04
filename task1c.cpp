// task1.cpp
//

/*Task 1: Sum of Random Variables. For Sn = L1 +    + Ln, assume Li are
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

using namespace std;

double Andersen_sidenius_Basu(int k, int n) {
	//double P_S_j_K = 0, P_S_j_k;
	if ((k > n) | (n < 0) | (k < 0)) {
		cout << "0 case : n=" << n << " and k = " << k << endl;
		return 0;
	}
	else {
		if ((n == 0) & (k == 0)) {
			//cout << "1 case" << endl;
			return 1;
		}
		else {
			if (k == n) {
				double Pk = 1;
				for (unsigned int i = 0; i <= n; i++) {
					double P_L = 1 - ((double)n / ((double)n + 1));
					Pk *= P_L;
				}
				return Pk;
			}
			else {
				double P_L = 1 - ((double)n / ((double)n + 1));
				cout << P_L << "k = " << k << " ; n =" << n << endl;
				return Andersen_sidenius_Basu(k, n - 1) * (1 - P_L) + Andersen_sidenius_Basu(k - 1, n - 1) * (P_L);
			}
		}
	}
}


double Hull_White(int n, int k) {
	double pi_n = 1;
	double ci;
	vector<double> Vn;
	vector<double> Un;
	for (unsigned int i = 0; i < k; i++) {
		cout << "k = " << i << endl;
		double Vi = 0;
		for (unsigned int j = 1; j <= n; j++) {
			double PL = 1 - ((double)j / ((double)n + 1));
			ci = PL / (1 - PL);
			cout << "P = " << PL << endl;
			cout << "ci = " << ci << endl;
			Vi += pow(ci, i);
		}
		Vn.push_back(Vi);
		double Ui = 0;
		if (Un.size() == 0) {
			Un.push_back(1);
		}
		
		for (unsigned int l = 0; l < Vn.size(); l++) {
			cout << "V" << l << " = " << Vn[l]<<endl;
			cout << "U" << Un.size() - l-1 << " = " << Un[Un.size() - l-1] << endl;
			Ui += pow(-1, l) * (Vn[l] * Un[Un.size() - l-1]);
		}
		Ui = Ui / (Un.size() + 1);
		Un.push_back(Ui);	
	}

	for (unsigned int i = 1; i <= n; i++) {
		double PL = 1 - ((double)i / ((double)n + 1));
		pi_n *= (1 - PL);
	}
	return pi_n*Un[Un.size()-1];
}


int main()
{
	///// Generate n independent variables Ui uniform on [0,1]
	// QQQ : Should we use congrential generator ? C++ tool ? 

	// With C++ tool
	default_random_engine generator;
	uniform_real_distribution<double> distribution(0.0, 1.0);

	// Generate Sn
	int n = 1000; // number of variables 

	double U_i;
	double p_i;
	double L_i;
	double S_n = 0;

	for (unsigned int i = 0; i < n; i++) {
		U_i = distribution(generator);
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
	cout << S_n << endl;

	// Compute Sn distribution
	int k = 500;
	double proba_k;

	/*for (unsigned int i = 0; i < n; i++) {

	}*/

	//  Andersen-Sidenius-Basu algorithm
	// Test
	cout << Andersen_sidenius_Basu(1, 1) << endl;


	// Hull-White algorithm
	cout << Hull_White(1, 1)<<endl;

	return 0;

}

// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.

// Exécuter le programme : Ctrl+F5 ou menu Déboguer > Exécuter sans débogage
// Déboguer le programme : F5 ou menu Déboguer > Démarrer le débogage

// Astuces pour bien démarrer : 
//   1. Utilisez la fenêtre Explorateur de solutions pour ajouter des fichiers et les gérer.
//   2. Utilisez la fenêtre Team Explorer pour vous connecter au contrôle de code source.
//   3. Utilisez la fenêtre Sortie pour voir la sortie de la génération et d'autres messages.
//   4. Utilisez la fenêtre Liste d'erreurs pour voir les erreurs.
//   5. Accédez à Projet > Ajouter un nouvel élément pour créer des fichiers de code, ou à Projet > Ajouter un élément existant pour ajouter des fichiers de code existants au projet.
//   6. Pour rouvrir ce projet plus tard, accédez à Fichier > Ouvrir > Projet et sélectionnez le fichier .sln.
