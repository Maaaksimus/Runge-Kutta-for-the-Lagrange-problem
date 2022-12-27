#include "mylib.h"

int main() {
	FILE *out, *func;
    vector<double> h;
	vector<vector<double> > x;
	double C[2], alph = ALPHA, t, s = 0;

	// shooting(C, alph);
    // cout << C[0] << " " << C[1] << endl;
	cout << "t1" << endl;
    // RungeKutta(x, C[0], C[1], alph);
    cout << "C at the start: " << C_1_default << " " << C_2_default << endl;
    RungeKutta(x, C_1_default, C_2_default, alph, h);

    // for (int i = 0; i < 4; i ++) {
    //     cout << x[x.size() - 1][i] << " ";
    // }

    for (int i = 0; i < x.size() - 1; i ++) {
        s += (0.5 * (x[i][0] + x[i + 1][0]) + pow(0.5 * (x[i][2] + x[i + 1][2]), 2)) * h[i];
    }

    cout << "in 0: " << x[0][1] << " " << x[0][3] << " " << endl;
    cout << "in pi/2: " << x[x.size() - 1][1] << " " << x[x.size() - 1][3] << " " << endl;
    cout << "Int: " << s << endl;
    // t = M_PI / 2.;
    // cout << endl << (-1/2.)*(pow(t,4) / 24. - M_PI / 4. * pow(t,3) / 6. + C_1_default*(-2)*t) << endl;

	return 0;
}

MyVector f(MyVector k, double alph) {
	MyVector KK;
    KK.x[0] = k.x[1];
	KK.x[1] = k.x[2];
	KK.x[2] = (-1)*k.x[3];
	KK.x[3] = exp(-alph * pow(k.x[0], 2)) * (1 - 2 * alph * pow(k.x[0], 2)) / 2.;
    // KK.print();
    return KK;
}

void countK(MyVector *k, double h, MyVector y, double alph) {  
    k[0] = h*f(y, alph);
    k[1] = h*f(y + (1 / 2.) * k[0], alph);
    k[2] = h*f(y + (1 / 4.) * (k[0] + k[1]), alph);
    k[3] = h*f(y + (-1)*k[1] + 2*k[2], alph);
    k[4] = h*f(y + (1 / 27.) * (7*k[0] + 10*k[1] + k[3]), alph);
    k[5] = h*f(y + (1 / 625.) * (28*k[0] + (-125)*k[1] + 546*k[2] + 54*k[3] + (-378)*k[4]), alph);
}

MyVector vec2Vec(vector<double> v) {
    MyVector V(v[0], v[1], v[2], v[3]);
    return V;
}

void RungeKutta(vector<vector<double> > &x, double C1, double C2, double alph, vector<double> &H) {
    ofstream out("h.txt");
    MyVector k_full[6], k_half_1[6], k_half_2[6],  buf(0, C1, 0, C2), E, k[6];
    vector<double> u;
    double pos = 0;
	double h = 1. / 100.;
    int step = 0;
    bool sw = 0, fin = 0;

    x.push_back(u);
	
	x[0].push_back(0);
    x[0].push_back(C1);
    x[0].push_back(0);
    x[0].push_back(C2);
    // cout << C1 << " " << C2 << endl;
    
    // while(M_PI / 2. - pos > EPS) {

    do {
        step ++;
        sw = 0;
        // cout << "t2" << endl;

        // cout << x.size() - 1 << endl;
        buf = vec2Vec(x[x.size() - 1]);

        // cout << "t3" << endl;

        do {

            if (pos + h > M_PI / 2.) {
                h = M_PI / 2. - pos;
                fin = 1;
            }

            countK(k_full, h, buf, alph);

            // k_full[0].print();
            
            // k[0] = h*f(pos, buf, alph);
            // k[1] = h*f(pos + 1 / 2. * h, buf + (1 / 2.) * k[0], alph);
            // k[2] = h*f(pos + 1 / 2. * h, buf + (1 / 4.) * (k[0] + k[1]), alph);
            // k[3] = h*f(pos + h, buf + (-1)*k[1] + 2*k[2], alph);
            // k[4] = h*f(pos + 2 / 3. * h, buf + (1 / 27.) * (7*k[0] + 10*k[1] + k[3]), alph);
            // k[5] = h*f(pos + 1 / 5. * h, buf + (1 / 625.) * (28*k[0] + (-125)*k[1] + 546*k[2] + 54*k[3] + (-378)*k[4]), alph);

            k[0] = h*f(buf, alph);
            k[1] = h*f(buf + (1 / 2.) * k[0], alph);
            k[2] = h*f(buf + (1 / 4.) * (k[0] + k[1]), alph);
            k[3] = h*f(buf + (-1)*k[1] + 2*k[2], alph);
            k[4] = h*f(buf + (1 / 27.) * (7*k[0] + 10*k[1] + k[3]), alph);
            k[5] = h*f(buf + (1 / 625.) * (28*k[0] + (-125)*k[1] + 546*k[2] + 54*k[3] + (-378)*k[4]), alph);

            // cout << "t4" << endl;
            // k[0].print();
            // k[1].print();
            // k[2].print();
            // k[3].print();
            // k[4].print();
            // k[5].print();


            // E = (1 / 336.) * ((-42)*k_full[0] + (-224)*k_full[2] + (-21)*k_full[3] + 162*k_full[4] + 125*k_full[5]);
            E = (1 / 336.) * ((-42)*k[0] + (-224)*k[2] + (-21)*k[3] + 162*k[4] + 125*k[5]);

            if (fin == 0) { 
                if (E.norm() < EPS / K) {
                    h *= 2;
                } else if (E.norm() > EPS) {
                    h /= 2;
                } else {
                    sw = 1;
                }
            } else {
                sw = 1;
            }

            out << "h: " << h << " " << E.x[0] << " "<< E.x[1] << " "<< E.x[2] << " "<< E.x[3] << " pos " << pos << endl;


        } while (sw == 0);

        // cout << "t6" << endl;

        vector<double> v;
        x.push_back(v);

        H.push_back(h);

        for (int i = 0; i < 4; i ++) {
            // x[step].push_back(x[step - 1][i] + 1. / 336. * (14*k_full[0].x[i] + 35*k_full[3].x[i] + 162*k_full[4].x[i] + 125*k_full[5].x[i]));
            x[step].push_back(x[step - 1][i] + 1. / 336. * (14*k[0].x[i] + 35*k[3].x[i] + 162*k[4].x[i] + 125*k[5].x[i]));
            // x[step].push_back(x[step - 1][i] + (1. / 6.)*(k[0].x[i] + 4*k[2].x[i] + k[3].x[i]));
        }

        // cout << "priv pos " << pos << endl;
        pos += h;
        // cout << "pos " << pos << endl;
    } while (fin == 0);
    // cout << "out" << endl;
}

void shooting(double *C, double alph) { // method "ready"
    
    vector<vector<double> > x_from_l;
    vector<vector<double> > x_from_r;
    vector<vector<double> > x;
    vector<double> h;
    double Jac_F[2][2], prev[2], a = 0, buf, det;

    C[0] = C_1_default;
    C[1] = C_2_default;

    for (int i = 1; i < 10; i ++) {

        a += alph / 9 * i;
        
        do {
            prev[0] = C[0];
            prev[1] = C[1];

            x.clear();
            x_from_l.clear();
            x_from_r.clear();
            h.clear();

            RungeKutta(x_from_r, C[0] - DIFF, C[1], a, h);
            // cout << "where are you?" << endl;
            RungeKutta(x_from_l, C[0] + DIFF, C[1], a, h);
            // cout << "here?" << endl;
            Jac_F[0][0] = (x_from_l[x_from_l.size() - 1][0] - x_from_r[x_from_r.size() - 1][0]) / DIFF / 2.;
            Jac_F[0][1] = (x_from_l[x_from_l.size() - 1][2] - x_from_r[x_from_r.size() - 1][2]) / DIFF / 2.;

            x_from_l.clear();
            x_from_r.clear();

            RungeKutta(x_from_r, C[0], C[1] - DIFF, a, h);
            RungeKutta(x_from_l, C[0], C[1] + DIFF, a, h);
            Jac_F[1][0] = (x_from_l[x_from_l.size() - 1][0] - x_from_r[x_from_r.size() - 1][0]) / DIFF / 2.;
            Jac_F[1][1] = (x_from_l[x_from_l.size() - 1][2] - x_from_r[x_from_r.size() - 1][2]) / DIFF / 2.;

            // invJac(Jac_F);

            det = Jac_F[0][0]*Jac_F[1][1] - Jac_F[0][1]*Jac_F[1][0];
            buf = Jac_F[0][0];

            Jac_F[0][0] = Jac_F[1][1] / det;
            Jac_F[0][1] *= -1 / det;
            Jac_F[1][0] *= -1 / det;
            Jac_F[1][1] = buf / det;

            RungeKutta(x, C[0], C[1], a, h);

            C[0] -= Jac_F[0][0]*(x[x.size() - 1][0] - 1) + Jac_F[0][1]*x[x.size() - 1][2];
            C[1] -= Jac_F[1][0]*(x[x.size() - 1][0] - 1) + Jac_F[1][1]*x[x.size() - 1][2];

            cout << C[0] << " " << C[1] << endl;
        } while(sqrt((C[0] - prev[0])*(C[0] - prev[0]) + (C[1] - prev[1])*(C[1] - prev[1])) > 1e-10);
    }
    x.clear();
    x_from_l.clear();
    x_from_r.clear();
}

void invJac(double J[2][2]) {
    double buf;
    // cout << "hi" << endl;
    // cout << J[0][0] << endl;
    // cout << "hello" << endl;
    buf = J[0][0];

    J[0][0] = J[1][1];
    J[0][1] *= -1;
    J[1][0] *= -1;
    J[1][1] = buf;
}