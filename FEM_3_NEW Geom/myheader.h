using namespace std;
using namespace Eigen;

tuple<double, double, double, MatrixXd, MatrixXd>get_Stress_Strains(double eps_p_m, double eps, double eta, double h, double t_total, double sigma_0, double E);