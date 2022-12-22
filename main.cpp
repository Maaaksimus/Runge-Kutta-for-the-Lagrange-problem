#include "mylib.h"
#include "MyVector.h"

int main() {
	FILE *out, *func;
	vector<vector<double> > x(4);
	double h, err_rate = 0, solve_rate = 0, C[2], alph = ALPHA;
	
	shooting(C, alph);
	RungeKutta(x, C[0], C[1], alph);

	return 0;
}

MyVector f(double t, MyVector k, double alph) {
	MyVector KK;
    KK.x[0] = k.x[1];
	KK.x[1] = k.x[2];
	KK.x[2] = -k.x[3];
	KK.x[3] = exp(-ALPHA * pow(k.x[0], 2)) * (1 - 2 * ALPHA * pow(k.x[0], 2));
    return KK;
}

MyVector vec2Vec(vector<double> v) {
    MyVector V;
    V.x[0] = v[0];
    V.x[1] = v[1];
    V.x[2] = v[2];
    V.x[3] = v[3];
    return V;
}

void RungeKutta(vector<vector<double> > &x, double C1, double C2, double alph) {
    MyVector k[6], buf(0, C1, 0, C2), E;
    double pos = 0;
	double h = 1. / 100.;
    int step = 0;
    bool sw = 0;
	
	x[0].push_back(0);
    x[0].push_back(C1);
    x[0].push_back(0);
    x[0].push_back(C2);
    
    while(M_PI / 2. - pos > EPS) {

        step ++;
        sw = 0;

        buf = vec2Vec(x[x.size() - 1]);

        do {
            k[0] = f(pos, buf, alph);
            k[1] = f(pos + 1 / 2. * h, buf + (1 / 2.) * k[0], alph);
            k[2] = f(pos + 1 / 2. * h, buf + (1 / 4.) * (k[0] + k[1]), alph);
            k[3] = f(pos + h, buf + (-1)*k[1] + 2*k[2], alph);
            k[4] = f(pos + 2 / 3. * h, buf + (1 / 27.) * (7*k[0] + 10*k[1] + k[3]), alph);
            k[5] = f(pos + 1 / 5. * h, buf + (1 / 625.) * (28*k[0] + -125*k[1] + 546*k[2] + 54*k[3] + -378*k[4]), alph);

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

            if (pos + h > M_PI / 2.) {
                h = M_PI / 2. - pos;
                sw = 0;
            }

        } while (sw == 0);

        vector<double> v;
        x.push_back(v);

        for (int i = 0; i < 4; i ++) {
            vector<double> v;
            x.push_back(v);
            x[step].push_back(x[step - 1][i] + 1. / 336. * (14*k[1].x[i] + 35*k[4].x[i] + 162*k[5].x[i] + 125*k[6].x[i]));
        }

        pos += h;
    }
}

void shooting(double *C, double alph) { // method "ready"
    
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

            RungeKutta(x_from_r, C[0] - EPS, C[1], a);
            RungeKutta(x_from_l, C[0] + EPS, C[1], a);
            Jac_F[0][0] = (x_from_l[0][x_from_l.size() - 1] - x_from_r[0][x_from_r.size() - 1]) / EPS / 2.;
            Jac_F[0][1] = (x_from_l[2][x_from_l.size() - 1] - x_from_r[2][x_from_r.size() - 1]) / EPS / 2.;

            for (int i = 0; i < 4; i ++) {
                x_from_l[i].clear();
                x_from_r[i].clear();
            }

            RungeKutta(x_from_r, C[0], C[1] - EPS, a);
            RungeKutta(x_from_l, C[0], C[1] + EPS, a);
            Jac_F[1][0] = (x_from_l[0][x_from_l.size() - 1] - x_from_r[0][x_from_r.size() - 1]) / EPS / 2.;
            Jac_F[1][1] = (x_from_l[2][x_from_l.size() - 1] - x_from_r[2][x_from_r.size() - 1]) / EPS / 2.;

            invJac(Jac_F);
            RungeKutta(x, C[0], C[1], a);

            C[0] -= Jac_F[0][0]*x[0][x.size() - 1] + Jac_F[0][1]*x[2][x.size() - 1];
            C[1] -= Jac_F[1][0]*x[0][x.size() - 1] + Jac_F[1][1]*x[2][x.size() - 1];
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