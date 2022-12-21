#include "mylib.h"

int main() {
	FILE *out, *func;
	vector<vector<double> > x(4);
	int n = 20, max_rate = 0;
	char *endptr;
	double h, err_rate = 0, solve_rate = 0, C[2], alph;

	printf("Enter the number of points: ");
	scanf("%d", &n);

	h = 1. / (n - 1);
	
	shooting(n, C, alph);
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

Vector f(double t, Vector k) {
	Vector KK;
    KK.x1 = k.x2;
	KK.x2 = k.x3;
	KK.x3 = -k.x4;
	KK.x4 = exp(-ALPHA * pow(k.x1, 2)) * (1 - 2 * ALPHA * pow(k.x1, 2));
    return KK;
}

Vector vec2Vec(vector<double> v) {
    Vector V;
    V.x1 = v[0];
    V.x2 = v[1];
    V.x3 = v[2];
    V.x4 = v[3];
    return V;
}

void RungeKutta(int n, vector<vector<double> > &x, double C1, double C2, double alph) {
    Vector k[6], buf(0, C1, 0, C2), E;
    double pos = 0;
	double h = 1. / 100.;
    int step = 0;
    bool sw = 0;
	
	x[0].push_back(0);
    x[1].push_back(C1);
    x[2].push_back(0);
    x[3].push_back(C2);
    
    while(pos + h < M_PI / 2.) { // надо заменить в классе x1, x2,... на x[i]

        step ++;

        buf = vec2Vec(x[x.size() - 1]);

        do {
            k[0] = f(pos, buf);
            k[1] = f(pos + 1 / 2. * h, buf + (1 / 2.) * k[0]);
            k[2] = f(pos + 1 / 2. * h, buf + (1 / 4.) * (k[0] + k[1]));
            k[3] = f(pos + h, buf + (-1)*k[1] + 2*k[2]);
            k[4] = f(pos + 2 / 3. * h, buf + (1 / 27.) * (7*k[0] + 10*k[1] + k[3]));
            k[5] = f(pos + 1 / 5. * h, buf + (1 / 625.) * (28*k[0] + -125*k[1] + 546*k[2] + 54*k[3] + -378*k[4]));

            E = (1 / 336.) * (-42*k[0] + -224*k[2] + -21*k[3] + 162*k[5] + 125*k[6]);

            for (int i = 0; i < 4; i ++) {
                if (E.x[i] < EPS / K) {
                    h *= 2;
                    break;
                } else if (E.x[i] > EPS) {
                    h /= 2;
                    break;
                } else {
                    sw = 1;
                }
            }
        } while (sw == 0);
        
    }
}

void shooting(int n, double *C, double alph) { // method "ready"
    
    vector<vector<double> > x_from_l(4);
    vector<vector<double> > x_from_r(4);
    vector<vector<double> > x(4);
    double Jac_F[2][2], prev[2], a = 0;

    C[0] = C_1_default;
    C[1] = C_2_default;

    for (int i = 1; i < 10; i ++) {

        a += alph / 9 * i;
        
        do {
            prev[0] = C[0];
            prev[1] = C[1];

            RungeKutta(n, x_from_r, C[0] - EPS, C[1], a);
            RungeKutta(n, x_from_l, C[0] + EPS, C[1], a);
            Jac_F[0][0] = (x_from_l[0][n - 1] - x_from_r[0][n - 1]) / EPS / 2.;
            Jac_F[0][1] = (x_from_l[2][n - 1] - x_from_r[2][n - 1]) / EPS / 2.;

            for (int i = 0; i < 4; i ++) {
                x_from_l[i].clear();
                x_from_r[i].clear();
            }

            RungeKutta(n, x_from_r, C[0], C[1] - EPS, a);
            RungeKutta(n, x_from_l, C[0], C[1] + EPS, a);
            Jac_F[1][0] = (x_from_l[0][n - 1] - x_from_r[0][n - 1]) / EPS / 2.;
            Jac_F[1][1] = (x_from_l[2][n - 1] - x_from_r[2][n - 1]) / EPS / 2.;

            invJac(Jac_F);
            RungeKutta(n, x, C[0], C[1], a);

            C[0] -= Jac_F[0][0]*x[0][n - 1] + Jac_F[0][1]*x[2][n - 1];
            C[1] -= Jac_F[1][0]*x[0][n - 1] + Jac_F[1][1]*x[2][n - 1];
        } while(sqrt((C[0] - prev[0])*(C[0] - prev[0]) + (C[1] - prev[1])*(C[1] - prev[1])) > 1e-10);
    }
}

void invJac(double J[2][2]) {
    double buf;

    buf = J[0][0];

    J[0][0] = J[1][1];
    J[0][1] *= -1;
    J[1][0] *= -1;
    J[1][1] = buf;
}