#include "mylib.h"

int main() {
    ofstream x1("x1.txt"), x3("x3.txt"), t_x1("targ1.txt"), t_x3("targ3.txt");
    vector<double> h;
	vector<vector<double> > x;
	double C[2], alph = 0, t, s = 0, pos = 0, x_1, x_3, s2 = 0;

	// shooting(C, alph);
    // cout << C[0] << " " << C[1] << endl;
	cout << "t1" << endl;
    // RungeKutta(x, C[0], C[1], alph, h);
    cout << "C at the start: " << C_1_default << " " << C_2_default << endl;
    RungeKutta(x, C_1_default, C_2_default, alph, h);

    // for (int i = 0; i < 4; i ++) {
    //     cout << x[x.size() - 1][i] << " ";
    // }

    for (int i = 0; i < x.size() - 1; i ++) {
        t = pos;
        s += (0.5 * (x[i][0] + x[i + 1][0]) + pow(0.5 * (x[i][2] + x[i + 1][2]), 2)) * h[i];
        
        // x_1 = (-1/2.)*(pow(t,4) / 24. - gig / 4. * pow(t,3) / 6. + C_1_default*(-2)*t);
        // x_3 = (-1/2.)*(pow(t,2) / 2. - gig / 4. * t);
        
        x1 << pos << " " << x[i][0] << endl;
        x3 << pos << " " << x[i][2] << endl;
        t_x1 << pos << " " << (-1/2.)*(pow(t,4) / 24. - gig / 4. * pow(t,3) / 6. + C_1_default*(-2)*t) << endl;
        t_x3 << pos << " " << (-1/2.)*(pow(t,2) / 2. - gig / 4. * t) << endl;
        pos += h[i];
    }

    cout << "in 0: " << x[0][1] << " " << x[0][3] << " " << endl;
    cout << "in pi/2: " << x[x.size() - 1][1] << " " << x[x.size() - 1][3] << " " << endl;
    cout << "Int: " << s << endl;
    // t = gig / 2.;
    // cout << (-1/2.)*(pow(t,3) / 6. - gig / 4. * pow(t,2) / 2. + C_1_default*(-2)) << endl;

	return 0;
}

MyVector f(MyVector k, double alph) {
	MyVector KK;
    KK.x[0] = k.x[1];
	KK.x[1] = k.x[2];
	KK.x[2] = (-1)*k.x[3];
	KK.x[3] = exp(-alph * pow(k.x[0], 2)) * (1 - 2 * alph * pow(k.x[0], 2)) / 2.;
    // KK.x[3] = 1 / 2.;
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
    MyVector k_full[6], k_half_1[6], k_half_2[6],  buf(0, C1, 0, C2), E, k[6], full, half, buf2;
    vector<double> v;
    double pos = 0;
	double h = 1. / 100000.;
    int step = 0;
    bool sw = 0, fin = 0;

    x.push_back(v);
	
	x[0].push_back(0)   ;
    x[0].push_back(C1);
    x[0].push_back(0);
    x[0].push_back(C2);
    // cout << C1 << " " << C2 << endl;
    
    // while(gig / 2. - pos > EPS) {

    do {
        step ++;
        sw = 0;
        // cout << "t2" << endl;

        // cout << x.size() - 1 << endl;
        buf = vec2Vec(x[x.size() - 1]);
        buf.print();

        // cout << "t3" << endl;

        do {

            if (pos + h > gig / 2.) {
                h = gig / 2. - pos;
                fin = 1;
            }

            // countK(k_full, h, buf, alph);

            // k_full[0].print();
            
            // k[0] = h*f(pos, buf, alph);
            // k[1] = h*f(pos + 1 / 2. * h, buf + (1 / 2.) * k[0], alph);
            // k[2] = h*f(pos + 1 / 2. * h, buf + (1 / 4.) * (k[0] + k[1]), alph);
            // k[3] = h*f(pos + h, buf + (-1)*k[1] + 2*k[2], alph);
            // k[4] = h*f(pos + 2 / 3. * h, buf + (1 / 27.) * (7*k[0] + 10*k[1] + k[3]), alph);
            // k[5] = h*f(pos + 1 / 5. * h, buf + (1 / 625.) * (28*k[0] + (-125)*k[1] + 546*k[2] + 54*k[3] + (-378)*k[4]), alph);

            k_full[0] = h*f(buf, alph);
            k_full[1] = h*f(buf + (1 / 2.) * k_full[0], alph);
            k_full[2] = h*f(buf + (1 / 4.) * (k_full[0] + k_full[1]), alph);
            k_full[3] = h*f(buf + (-1)*k_full[1] + 2*k_full[2], alph);
            k_full[4] = h*f(buf + (1 / 27.) * (7*k_full[0] + 10*k_full[1] + k_full[3]), alph);
            k_full[5] = h*f(buf + (1 / 625.) * (28*k_full[0] + (-125)*k_full[1] + 546*k_full[2] + 54*k_full[3] + (-378)*k_full[4]), alph);

            k_full[0] = h*f(buf, alph);
            k_full[1] = h*f(buf + (1 / 2.) * k_full[0], alph);
            k_full[2] = h*f(buf + (1 / 4.) * (k_full[0] + k_full[1]), alph);
            k_full[3] = h*f(buf + (-1)*k_full[1] + 2*k_full[2], alph);
            k_full[4] = h*f(buf + (1 / 27.) * (7*k_full[0] + 10*k_full[1] + k_full[3]), alph);
            k_full[5] = h*f(buf + (1 / 625.) * (28*k_full[0] + (-125)*k_full[1] + 546*k_full[2] + 54*k_full[3] + (-378)*k_full[4]), alph);

/*
            k_half_1[0] = (h / 2.)*f(buf, alph);
            k_half_1[1] = (h / 2.)*f(buf + (1 / 2.) * k_half_1[0], alph);
            k_half_1[2] = (h / 2.)*f(buf + (1 / 4.) * (k_half_1[0] + k_half_1[1]), alph);
            k_half_1[3] = (h / 2.)*f(buf + (-1)*k_half_1[1] + 2*k_half_1[2], alph);
            k_half_1[4] = (h / 2.)*f(buf + (1 / 27.) * (7*k_half_1[0] + 10*k_half_1[1] + k_half_1[3]), alph);
            k_half_1[5] = (h / 2.)*f(buf + (1 / 625.) * (28*k_half_1[0] + (-125)*k_half_1[1] + 546*k_half_1[2] + 54*k_half_1[3] + (-378)*k_half_1[4]), alph);

            buf2 = buf + 1. / 336. * (14*k_half_1[0] + 35*k_half_1[3] + 162*k_half_1[4] + 125*k_half_1[5]);

            k_half_2[0] = (h / 2.)*f(buf2, alph);
            k_half_2[1] = (h / 2.)*f(buf2 + (1 / 2.) * k_half_2[0], alph);
            k_half_2[2] = (h / 2.)*f(buf2 + (1 / 4.) * (k_half_2[0] + k_half_2[1]), alph);
            k_half_2[3] = (h / 2.)*f(buf2 + (-1)*k_half_2[1] + 2*k_half_2[2], alph);
            k_half_2[4] = (h / 2.)*f(buf2 + (1 / 27.) * (7*k_half_2[0] + 10*k_half_2[1] + k_half_2[3]), alph);
            k_half_2[5] = (h / 2.)*f(buf2 + (1 / 625.) * (28*k_half_2[0] + (-125)*k_half_2[1] + 546*k_half_2[2] + 54*k_half_2[3] + (-378)*k_half_2[4]), alph);
            
            // cout << "t4" << endl;
            // k[0].print();
            // k[1].print();
            // k[2].print();
            // k[3].print();
            // k[4].print();
            // k[5].print();

            full = buf + 1. / 336. * (14*k_full[0] + 35*k_full[3] + 162*k_full[4] + 125*k_full[5]);
            half= buf2 + 1. / 336. * (14*k_half_2[0] + 35*k_half_2[3] + 162*k_half_2[4] + 125*k_half_2[5]);

*/
            E = (1 / 336.) * ((-42)*k_full[0] + (-224)*k_full[2] + (-21)*k_full[3] + 162*k_full[4] + 125*k_full[5]);
/*
            if (fin == 0) { 
                if ((full + (-1)*half).norm() < EPS / K) {
                // if (E.norm() < EPS / K) {
                    h *= 2;
                } else if ((full + (-1)*half).norm() > EPS) {
                // } else if (E.norm() > EPS) {
                    h /= 2;
                } else {
                    sw = 1;
                }
            } else {
                sw = 1;
            }
            */

            out << "h: " << h << " " << E.x[0] << " "<< E.x[1] << " "<< E.x[2] << " "<< E.x[3] << " pos " << pos << " step: " << step << endl;
            sw = 1;

        } while (sw == 0);

        // cout << "t6" << endl;

        x.push_back(v);

        H.push_back(h);

        buf2 = buf + (1./24.)*k_full[0] + (5./48.)*k_full[3] + (27./56.)*k_full[4] + (125./336.)*k_full[5];

        for (int i = 0; i < 4; i ++) {
            // x[step].push_back(x[step - 1][i] + 1. / 336. * (14*k_full[0].x[i] + 35*k_full[3].x[i] + 162*k_full[4].x[i] + 125*k_full[5].x[i]));
            // x[step].push_back(x[step - 1][i] + 1. / 336. * (14*k_full[0].x[i] + 35*k_full[3].x[i] + 162*k_full[4].x[i] + 125*k_full[5].x[i]));
            // x[step].push_back(x[step - 1][i] + (1./24.)*k_full[0].x[i] + (5./48.)*k_full[3].x[i] + (27./56.)*k_full[4].x[i] + (125./336.)*k_full[5].x[i]);
            // x[step].push_back(x[step - 1][i] + (1. / 6.)*(k[0].x[i] + 4*k[2].x[i] + k[3].x[i]));
            x[step].push_back(buf2.x[i]);
        }

        // cout << "priv pos " << pos << endl;
        pos += h;
        // cout << "pos " << pos << endl;
    } while (fin == 0);
    // printf("pos = %.10f\n", pos);//10.677
    // printf("pi = %.20lf\n", gig);//10.677
    // cout << "with changes" << endl;
    // cout << setprecision(20) << gig << "pa pam" << endl;;
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

        a = alph / 9. * i;
        cout << "a: " << a << endl;
        
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
            Jac_F[1][0] = (x_from_l[x_from_l.size() - 1][2] - x_from_r[x_from_r.size() - 1][2]) / DIFF / 2.;

            x_from_l.clear();
            x_from_r.clear();

            RungeKutta(x_from_r, C[0], C[1] - DIFF, a, h);
            RungeKutta(x_from_l, C[0], C[1] + DIFF, a, h);
            Jac_F[0][1] = (x_from_l[x_from_l.size() - 1][0] - x_from_r[x_from_r.size() - 1][0]) / DIFF / 2.;
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
        } while(sqrt((C[0] - prev[0])*(C[0] - prev[0]) + (C[1] - prev[1])*(C[1] - prev[1])) > 1e-4);
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