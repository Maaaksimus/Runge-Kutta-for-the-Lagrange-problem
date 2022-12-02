#include "mylib.h"

int main() {
	FILE *out, *func;
	vector<double> u;
	int n = 20, max_rate = 0;
	char *endptr;
	double h, err_rate = 0, solve_rate = 0, C[2];

	printf("Enter the number of points: ");
	scanf("%d", &n);

	h = 1. / (n - 1);
	
	// shooting(n, e);
	RungeKutta();

	out = fopen("res.txt", "w");
	func = fopen("targ.txt", "w");

	max_rate = 0;

	for (int i = 0; i < n; i ++) {
		fprintf(out, "%lf %lf\n", i*h, u[0][i]);
	}

	fclose(out);
	fclose(func);

	printf("Error rate: %le\n", sqrt(err_rate));
	printf("Check results in the file\n");
	return 0;
}

// double diff(double f(double), double x) {
//     return (f(x + EPS) - f(x)) / EPS;
// }

void RungeKutta(int n, vector<vector<double> > *x, double C1, double C2, double alph) {

}
void Newton(double *C, double alph) { // method "ready"
    
    vector<vector<double> > x_from_l(4);
    vector<vector<double> > x_from_r(4);
    vector<vector<double> > x(4);
    double Jac_F[2][2], prev[2];
    int n;

    C[0] = C_1_default;
    C[1] = C_2_default;

    do {
        prev[0] = C[0];
        prev[1] = C[1];

        RungeKutta(n, &x_from_r, C[0] - EPS, C[1], alph);
        RungeKutta(n, &x_from_l, C[0] + EPS, C[1], alph);
        Jac_F[0][0] = (x_from_l[0][n - 1] - x_from_r[0][n - 1]) / EPS;
        Jac_F[0][1] = (x_from_l[2][n - 1] - x_from_r[2][n - 1]) / EPS;

        for (int i = 0; i < 4; i ++) {
            x_from_l[i].clear();
            x_from_r[i].clear();
        }

        RungeKutta(n, &x_from_r, C[0], C[1] - EPS, alph);
        RungeKutta(n, &x_from_l, C[0], C[1] + EPS, alph);
        Jac_F[1][0] = (x_from_l[0][n - 1] - x_from_r[0][n - 1]) / EPS;
        Jac_F[1][1] = (x_from_l[2][n - 1] - x_from_r[2][n - 1]) / EPS;

        invJac(Jac_F);
        RungeKutta(n, &x, C[0], C[1], alph);

        C[0] -= Jac_F[0][0]*x[0][n - 1] + Jac_F[0][1]*x[2][n - 1];
        C[1] -= Jac_F[1][0]*x[0][n - 1] + Jac_F[1][1]*x[2][n - 1];
    } while(sqrt((C[0] - prev[0])*(C[0] - prev[0]) + (C[1] - prev[1])*(C[1] - prev[1])) > 1e-10);
}

void invJac(double J[2][2]) {
    double buf;

    buf = J[0][0];

    J[0][0] = J[1][1];
    J[0][1] *= -1;
    J[1][0] *= -1;
    J[1][1] = buf;
}