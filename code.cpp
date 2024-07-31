#include <iostream>
#include <vector>
#include <cmath>

using namespace std;


long double normal(long double x, long double mean, long double v) {
    long double a = -0.5 * pow((x - mean) / sqrt(v), 2);
    long double b = sqrt(2 * M_PI * v);
    long double c = exp(a) / b;
    return c;
}



void solve(const vector<long double>& returns, long double a00, long double a01, long double a10, long double a11, long double pi0, long double pi1) {
    int T = returns.size();
    vector<long double> mu = {0, 0};
    vector<long double> v = {1, 1}; 
    vector<vector<long double>> Al(T, vector<long double>(2)); 
    vector<vector<long double>> Be(T, vector<long double>(2)); 
    vector<vector<long double>> Y(T, vector<long double>(2));
    vector<vector<vector<long double>>> xi(T - 1, vector<vector<long double>>(2, vector<long double>(2))); 

    long double oa00, oa01, oa10, oa11, opi0, opi1;
    long double da00, da01, da10, da11, dpi0, dpi1,dmu0,dmu1,dv0,dv1;
    long double omu0=mu[0];
    long double omu1=mu[1];
    long double ov0=v[0];
    long double ov1=v[1];

    do {
            oa00 = a00;
            oa01 = a01; 
            oa10 = a10; 
            oa11 = a11;
            opi0 = pi0; 
            opi1 = pi1;
            Al[0][0] = pi0 * normal(returns[0], mu[0], v[0]);
            Al[0][1] = pi1 * normal(returns[0], mu[1], v[1]);
            Be[T - 1][0] = 1;
            Be[T - 1][1] = 1;
        

        for (int t = T - 2; t >= 0; --t) {
            Be[t][0] = a00 * normal(returns[t + 1], mu[0], v[0]) * Be[t + 1][0] + a01 * normal(returns[t + 1], mu[1], v[1]) * Be[t + 1][1];
            Be[t][1] = a10 * normal(returns[t + 1], mu[0], v[0]) * Be[t + 1][0] + a11 * normal(returns[t + 1], mu[1], v[1]) * Be[t + 1][1];
        }

        for (int t = 1; t < T; ++t) {
            Al[t][0] = Al[t - 1][0] * a00 * normal(returns[t], mu[0], v[0]) + Al[t - 1][1] * a10 * normal(returns[t], mu[0], v[0]);
            Al[t][1] = Al[t - 1][0] * a01 * normal(returns[t], mu[1], v[1]) + Al[t - 1][1] * a11 * normal(returns[t], mu[1], v[1]);
        }

        

        for (int t = 0; t < T; ++t) {
            long double sum_Y = 0;
            for (int i = 0; i < 2; ++i) {
                Y[t][i] = Al[t][i] * Be[t][i];
                sum_Y += Y[t][i];
            }
            for (int i = 0; i < 2; ++i) {
                Y[t][i] /= sum_Y;
            }
            if (t < T - 1) {
                long double sum_xi = 0;
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        xi[t][i][j] = Al[t][i] * (i == 0 ? a00 : a10) * Be[t + 1][j] * normal(returns[t + 1], mu[j], v[j]);
                        sum_xi += xi[t][i][j];
                    }
                }
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        xi[t][i][j] /= sum_xi;
                    }
                }
            }
            
        }
        
        long double na00 = 0, na01 = 0, na10 = 0, na11 = 0;
        
        long double sum_Y_0 = 0, sum_Y_1 = 0;
        long double sum_Y_y_0 = 0, sum_Y_y_1 = 0;

        for (int t = 0; t < T - 1; ++t) {
            na00 += xi[t][0][0];
            na01 += xi[t][0][1];
            na10 += xi[t][1][0];
            na11 += xi[t][1][1];
            sum_Y_0 += Y[t][0];
            sum_Y_1 += Y[t][1];
            sum_Y_y_0 += Y[t][0] * returns[t];
            sum_Y_y_1 += Y[t][1] * returns[t];
        }

        na00 /= sum_Y_0;
        na01 /= sum_Y_0;
        na10 /= sum_Y_1;
        na11 /= sum_Y_1;
        omu0=mu[0];
        omu1=mu[1];
        ov0=v[0];
        ov1=v[1];
        mu[0] = sum_Y_y_0 / sum_Y_0;
        mu[1] = sum_Y_y_1 / sum_Y_1;
        v[0]=0;
        v[1]=0;
        for (int t = 0; t < T; ++t) {
            
            v[0] += Y[t][0] * pow(returns[t] - mu[0], 2);
            v[1] += Y[t][1] * pow(returns[t] - mu[1], 2);
        }
        v[0] /= sum_Y_0;
        v[1] /= sum_Y_1;
        a00=na00;
        a01=na01;
        a10=na10;
        a11=na11;
        pi0 = Y[0][0];
        pi1 = Y[0][1];
        da00 = abs(oa00 - a00);
        da01 = abs(oa01 - a01);
        da10 = abs(oa10 - a10);
        da11 = abs(oa11 - a11);
        dpi0 = abs(opi0 - pi0);
        dpi1 = abs(opi1 - pi1);
        dv1=abs(ov1-v[1]);
        dv0=abs(ov0-v[0]);  
        dmu0=abs(omu0-mu[0]);
        dmu1=abs(omu1-mu[1]);
    } while (da00> 1e-8 || da01 > 1e-8 || da10 > 1e-8 ||
             da11 > 1e-8 || dpi0 > 1e-8 || dpi1 > 1e-8 || dmu1 > 1e-8 || dmu0 > 1e-8 || dv1 > 1e-8 || dv0 > 1e-8);

    // cout<<v[0]<<" "<<v[1]<<endl;
    //da00+da01+da10+da11+dpi0+dpi1+dmu1+dmu0+dv1+dv0>1e-8
        long double npi0,npi1;
        long double bull;
        long double bear;
    for (int t = 0; t < T; ++t) {
        
        npi0 = Y[t][0];
        npi1 = Y[t][1];
        // cout<<Y[t][1]<<endl;
        if(mu[0]>=mu[1]){
            bull = npi0;
            bear = npi1;
        }
        else{
            bull = npi1;
            bear = npi0;
        }
        

        if (bull > bear) {
            cout << "Bull" << endl;
        } else {
            cout << "Bear" << endl;
        }

    
    }

}

int main() {

    long double a00, a01, a10, a11, pi0, pi1;
    cin >> a00 >> a01 >> a10 >> a11 >> pi0 >> pi1;
    int T;
    cin >> T;
    vector<long double> returns;
    returns.resize(T,0);
    for (int i = 0; i < T; ++i) {
        long double ai;
        cin>>ai;
        returns[i]=ai;
    }
    solve(returns, a00, a01, a10, a11, pi0, pi1);
    

}
