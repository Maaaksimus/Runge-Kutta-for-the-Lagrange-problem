#include "mylib.h"

int main() {
    ofstream x1("x1.txt"), x3("x3.txt"), t_x1("targ1.txt"), t_x3("targ3.txt");
    vector<double> h;
	vector<vector<double> > x;
    vector<MyVector> X;
	double C[2], alph = ALPHA, t, s = 0, pos = 0;

	shooting(C, alph);
    cout << C[0] << " " << C[1] << endl;
    // cout << "C at the start: " << C_1_default << " " << C_2_default << endl;
    // RungeKutta(X, C_1_default, C_2_default, alph, h);
    RungeKutta(X, C[0], C[1], alph, h);

    for (int i = 0; i < X.size() - 1; i ++) {
        // s += (0.5 * (x[i][0] + x[i + 1][0]) + pow(0.5 * (x[i][2] + x[i + 1][2]), 2)) * h[i];
        // x1 << pos << " " << x[i][0] << endl;
        // x3 << pos << " " << x[i][2] << endl;
        // t = pos;
        // t_x1 << pos << " " << (-1/2.)*(pow(t,4) / 24. - M_PI_4 * pow(t,3) / 6. + C_1_default*(-2)*t) << endl;
        // t_x3 << pos << " " << (-1/2.)*(pow(t,2) / 2. - M_PI_4 * t) << endl;
        // pos += h[i];
        s += L(X[i], X[i+1], alph) * h[i];
    }

    cout << "in 0: "; X[0].print();
    cout << "in pi/2: "; X[X.size() - 1].print();
    cout << "Int: " << s << endl;
    // t = M_PI / 2.;
    // cout << (-1/2.)*(pow(t,3) / 6. - M_PI / 4. * pow(t,2) / 2. + C_1_default*(-2)) << endl;

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

void shooting(double *C, double alph) { // method "ready"
    vector<MyVector> XL, XR, X_curr;
    vector<double> h;
    double Jac_F[2][2], prev[2], a = 0, buf, det;

    C[0] = C_1_default;
    C[1] = C_2_default;

    for (int i = 1; i < 10; i ++) {

        a = alph / 9. * i;
        cout << "a: " << a << endl;
        
        do {
            prev[0] = C[0];
            prev[1] = C[1];

            X_curr.clear();
            XL.clear();
            XR.clear();
            h.clear();

            RungeKutta(XR, C[0] - DIFF, C[1], a, h);
            // cout << "where are you?" << endl;
            RungeKutta(XL, C[0] + DIFF, C[1], a, h);
            // cout << "here?" << endl;
            Jac_F[0][0] = (XL[XL.size() - 1].x[0] - XR[XR.size() - 1].x[0]) / DIFF / 2.;
            Jac_F[1][0] = (XL[XL.size() - 1].x[2] - XR[XR.size() - 1].x[2]) / DIFF / 2.;

            XL.clear();
            XR.clear();

            RungeKutta(XR, C[0], C[1] - DIFF, a, h);
            RungeKutta(XL, C[0], C[1] + DIFF, a, h);
            Jac_F[0][1] = (XL[XL.size() - 1].x[0] - XR[XR.size() - 1].x[0]) / DIFF / 2.;
            Jac_F[1][1] = (XL[XL.size() - 1].x[2] - XR[XR.size() - 1].x[2]) / DIFF / 2.;

            // invJac(Jac_F);

            det = Jac_F[0][0]*Jac_F[1][1] - Jac_F[0][1]*Jac_F[1][0];
            buf = Jac_F[0][0];

            Jac_F[0][0] = Jac_F[1][1] / det;
            Jac_F[0][1] *= -1 / det;
            Jac_F[1][0] *= -1 / det;
            Jac_F[1][1] = buf / det;

            RungeKutta(X_curr, C[0], C[1], a, h);

            C[0] -= Jac_F[0][0]*(X_curr[X_curr.size() - 1].x[0] - 1) + Jac_F[0][1]*X_curr[X_curr.size() - 1].x[2];
            C[1] -= Jac_F[1][0]*(X_curr[X_curr.size() - 1].x[0] - 1) + Jac_F[1][1]*X_curr[X_curr.size() - 1].x[2];

            // cout << C[0] << " " << C[1] << endl;
        } while(sqrt((C[0] - prev[0])*(C[0] - prev[0]) + (C[1] - prev[1])*(C[1] - prev[1])) > 1e-7);
    }

    X_curr.clear();
    XL.clear();
    XR.clear();
}

void RungeKutta(vector<MyVector> &X, double C1, double C2, double alph, vector<double> &H) {
    ofstream out("h.txt");
    MyVector k_full[6], k_half_1[6], k_half_2[6],  buf(0, C1, 0, C2), E, k[6], full, half, buf2;
    // vector<MyVector> X;
    double pos = 0;
	double h = 1. / 1000.;
    int step = 0;
    bool sw = 0, fin = 0, mult = 0;

    X.push_back(buf);

    do {
        step ++;
        sw = 0;

        do {

            if (pos + h > M_PI_2) {
                h = M_PI_2 - pos;
                fin = 1;
            }

            buf = X[step - 1];
          
            full = countStep(h, buf, alph);

            buf2 = countStep(h/2., buf, alph);
            half = countStep(h/2., buf2, alph);


            if (fin == 0) { 
                if ((full + (-1)*half).norm() < EPS / K) {
                // if (E.norm() < EPS / K) {
                    // cout << (full + (-1)*half).norm() << endl;
                    if (mult == 1) {break;}
                    h *= 2;
                    mult = 1;
                } else if ((full + (-1)*half).norm() > EPS) {
                // } else if (E.norm() > EPS) {
                    if (mult == 1) {h /= 2; break;}
                    h /= 2;
                    mult = 1;
                } else {
                    sw = 1;
                }
            } else {
                sw = 1;
            }

            // out << "h: " << h << " " << E.x[0] << " "<< E.x[1] << " "<< E.x[2] << " "<< E.x[3] << " pos " << pos << endl;
            out << "h: " << h << " pos: " << pos <<  " " << full.x[1] << " " << full.x[3] << " error: " << (full + (-1)*half).norm() << endl;
           
        } while (sw == 0);


        H.push_back(h);

        X.push_back(full);

        pos += h;

    } while (fin == 0);

    out.close();
}

MyVector countStep(double h, MyVector y, double alph) {
    MyVector k[6];

    k[0] = h*f(y, alph);
    k[1] = h*f(y + (1./2.)*(k[0]), alph);
    k[2] = h*f(y + (1./4.)*(k[0] + k[1]), alph);
    k[3] = h*f(y + (-1)*k[1] + 2*k[2], alph);
    k[4] = h*f(y + (1./27.)*(7*k[0] + 10*k[1] + k[3]), alph);
    k[5] = h*f(y + (1./625.)*(28*k[0] + (-125)*k[1] + 546*k[2] + 54*k[3] + (-378)*k[4]), alph);

    return y + (1./24.)*k[0] + (5./48.)*k[3] + (27./56.)*k[4] + (125./336.)*k[5];
}

double L(MyVector x_curr, MyVector x_next, double alph) {
    double x = (1./2.)*(x_curr.x[0] + x_next.x[0]);
    double x_pp = (1./2.)*(x_curr.x[2] + x_next.x[2]);
    
    return x_pp*x_pp + x*exp((-1)*alph*x*x);
}
