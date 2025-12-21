#include <iostream>
using namespace std;

// -- Class/Struct defenitions -- 

struct OptionParameters{

    long double S, K, T, r, sig;


    OptionParameters(long double S_, long double K_, long double T_, long double r_, long double sig_) {
        S = S_;     // Stock price
        K = K_;     // Strike price
        T = T_;     // Time until expiry
        r = r_;     // Risk free Rate of Interest
        sig = sig_; // Volatility
    }
    
    void print() const {
        cout << "S: " << S << ", K: " << K
            << ", T: " << T << ", r: " << r
            << ", sig: " << sig << endl;
    }
};
