#include<iostream>
#include<Eigen/Dense>
#include<iomanip>
#include<fstream>
#include"myheader.h"

using namespace std;
using namespace Eigen;

tuple<double, double, double, MatrixXd, MatrixXd>get_Stress_Strains(double eps_p_m, double eps, double eta, double h, double t_total, double sigma_0, double E) {
	int s = (t_total / h);
	MatrixXd sigmas(s, 0);
	MatrixXd ep_p_vals(s, 0);
	double Ct = E;
	double sigma_m, eps_p_next, eps_p_prev, sigma_next;

	for (int i = 0; i < s; i++) {
		sigma_m = E * (eps - eps_p_m);
		eps_p_next = eps_p_m + h * eta * signum(sigma_m) * Lambda(sigma_m, sigma_0);

		for (int j = 0; j < 1000; j++) {
			eps_p_prev = eps_p_next;
			sigma_next = E * (eps - eps_p_next);
			eps_p_next = eps_p_m + h * eta * signum(sigma_next) * Lambda(sigma_next, sigma_0);

			if (abs(eps_p_next - eps_p_prev) < 1e-5) break;
		}

		sigma_m = E * (eps - eps_p_next);
		eps_p_m = eps_p_next;

		sigmas(i, 0) = sigma_m;
		ep_p_vals(i, 0) = eps_p_next;
	}
	if (sigma_m > sigma_0) {
		Ct = (E * sigma_0) / (sigma_0 + E * eta * t_total);
	}
	else Ct = E;

	return { sigma_m, eps_p_next, Ct,sigmas, ep_p_vals };
}
