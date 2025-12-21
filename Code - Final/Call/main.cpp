
#include "OptionParameters.h" // made a header file to store classes
#include<iostream>
#include <fstream>
#include <cmath>

using namespace std;

struct OptionParameters;     
long double factorial(long double n);
long double COSCer(long double N, long double j, long double quu, long double qdd);
long double zScore(long double val);
long double BSM_Call(OptionParameters t);
long double BM_Call(OptionParameters params, long double StepDivisions = 1000.0);

//
//
/* 

The main() function:
We are trying to Create the main function to iterate through various values of a 
Hyper parameter (say Stock price (S)) and See how the parameter changes affect the Convergence of the model

*/
//
//

int main(){
    

    OptionParameters Params1(
        /*S*/   100.0,
        /*K*/   120.0,
        /*T*/   5.0,
        /*r*/   0.1,
        /*sig*/ 0.7
    );

    // save to csv for plotting purposes
    string base_file, file_name;
    base_file = "/Users/rrvsants/Documents/textfiles/qmul/Term1/PLA/0-final/";
    file_name = "-sig.csv"; //  "S.csv" "K.csv" "r.csv" "T.csv" "sig.csv"
    ofstream out(base_file + file_name);
    out << "Step,hp,BSM,Binomial,Spread\n";

    double StepInput, limiter = 1000.0, hpmax = 400.0;


    for(double hp=1.0; hp <= hpmax; hp++){
        // Params1.T = hp;
        // Params1.S = hp;
        // Params1.r = 0.99 * Params1.sig * hp / hpmax ;  // r / N < sig * srt(del t), is the case and so Su > So(1+R) > Sd
        // Params1.K = hp;
        Params1.sig = Params1.r * 1.2 + hp / hpmax * 1 ;
        
        long double CallPrice_BSM = BSM_Call(Params1);

        for(StepInput = 1.0; StepInput <= limiter; StepInput++){
            // cout<<"\n-------------------------------------------"<<endl;
            
            // Params1.S = StepInput;
            long double CallPrice_BM = BM_Call(Params1, StepInput);
            if(isnan(CallPrice_BM)){
                cout << "* * * * NaN detected at Step = " << StepInput << ". * * * *.\n";
                break;

            }

            long double diff = abs(CallPrice_BSM - CallPrice_BM);
            
            out << StepInput << ","
                << hp   << ","
                << CallPrice_BSM  << ","
                << CallPrice_BM << ","
                << diff << "\n";

            cout<<StepInput<<". "<<"\nBSM :"<<CallPrice_BSM<<"\n"
                <<"BM  :"<<CallPrice_BM<<endl<<"spread : "<<CallPrice_BSM - CallPrice_BM;

        }
    }

    out.close();
    return 0;
    
} 

// -- Function defenitions -- 

//
//
/* 

The sub functions:

BM_Call: this function calculates the the Call Option price using the N-period Binomial model setup, the input parameters are : S, K, T, sig, r, N and output => Ct
        Uses subfunction COSCer

BSM_Call: this function calculates the the Call Option price using the Black-Scholes Formula, the input parameters are : S, K, T, sig, r and output => Ct
        Uses the subfuncion zScore

    Sub Sub Functions:

    COSCer: This is a code to calculate each term of the binomial sample, so (N C j) * (u ^ j) * (d ^ (N - j)),
        but we know there is a limit to the range of numbers that c++ can store, so we use an optimal iteration to scale the 
        value back down everytime before the factorial values jump, so the name CRR Optimised Scaled Combinations
            > CRR - The function needs to output a part of the CRR expression that outputs the chance of the Step happening as an Binomial trial
            > Optimised Scaled - But, lucky for us u and d are < 1, so we try scaling down the values of each term in 1*2*...*j of j!
            > Combinations - We finally do it for each factorial in the (N C j) combination and arrive at the used expression, which allows us to go until 930 steps!

    zScore: We require the cumulative normal distrbution function that does not directly exist in the c++ libraries, so we use the error function to get to the cdf.

*/
//
//

long double BM_Call(OptionParameters params, long double StepDivisions){
    // Param.S
    // Param.K
    // Param.T
    // Param.r
    // Param.sig
    // qu, qd, u, d, qu_, qd_, j_min, StepDivisions
    // StepDivisions = N,while Param.T = Time periods of r executable

    // Continous per T period compunding -> discrete N Step compounding
    long double r_discrete = (exp(params.r * params.T / StepDivisions) - 1) * StepDivisions;

    long double R_bin = (r_discrete) / StepDivisions;
    
    // step size, and Risk Neutral Probabilities
    long double u, d, qu, qd, qu_, qd_;
    u = exp(params.sig * sqrt( params.T/ StepDivisions));
    d = exp(-params.sig * sqrt( params.T/ StepDivisions));

    qu = (1.0 + (R_bin) - d) / (u - d);
    qd = 1 - qu;

    // have not used them
    qu_ = (qu * u) / (1 + ((R_bin)));
    qd_ = 1 - qu_;

    long double j_min = (log(params.K) - log(params.S) - StepDivisions * log(d)) / log(u/d);
    

    //j_min checked
    j_min = j_min - fmod(j_min, 1.0) + 1.0;
    
    // check derived variables
    // cout <<"-------------------------------------------"<<endl
    //     << "u   = " << u   << endl
    //     << "d   = " << d   << endl
    //     << "j_min   = " << j_min   << endl
    //     << "Rate of Interest per step - MultiBinMod  = " << ((R_bin))   << endl
    //     << "Step Size: " << param.T / StepDivisions << endl
    //     << "qu  = " << qu  << endl
    //     << "qd  = " << qd  << endl
    //     << "qu_ = " << qu_ << endl
    //     << "qd_ = " << qd_ << endl
    //     <<"-------------------------------------------"<<endl;
    
    long double Ct, S_coeff = 0.0, K_coeff = 0.0;

    for(long double j = j_min; j <= StepDivisions ; j += 1.0){
        
        S_coeff += COSCer(StepDivisions, j, qu * u, qd * d);// *  pow(qu * u, j) * pow(qd * d, StepDivisions - j);

        K_coeff += COSCer(StepDivisions, j, qu, qd);// * pow(qu, j) * pow(qd, StepDivisions - j);
    }
    
    Ct = (params.S / pow(1 + R_bin, StepDivisions)) * S_coeff - (params.K / pow(1 + R_bin, StepDivisions)) * K_coeff;
    
    // check derived variables
    // cout <<"-------------------------------------------"<<endl
    //     << "\nT1   = " << S_coeff   << endl
    //     << "T2   = " << K_coeff   << endl
    //     << "S pv   = " << (params.S / pow(1 + R_bin, StepDivisions))   << endl
    //     << "K pv   = " << (params.K / pow(1 + R_bin, StepDivisions))   << endl
    //     << "S term   = " << (params.S / pow(1 + R_bin, StepDivisions)) * S_coeff   << endl
    //     << "K term   = " << (params.K / pow(1 + R_bin, StepDivisions)) * K_coeff   << endl
    //     << "Ct   = " << Ct   << endl
    //     <<"-------------------------------------------"<<endl;
    // cout<<"* * * * * * * * * * * * * * * * * * * * * * * * * * * * \n";



    return Ct;

}

// BSM: S, K, T, r ==> Ct
long double BSM_Call(OptionParameters params){ //long double S, long double K, long double T, long double r, long double sig
    long double d1, d2;

    // check derived variables
    
    
    d1 = (log(params.S / params.K) + (params.r + (params.sig*params.sig/2)) * params.T) / (params.sig * sqrt(params.T));
    d2 = d1 - params.sig * sqrt(params.T);
    
    long double Ct = params.S * zScore(d1) - params.K * zScore(d2) * exp(-params.r * params.T);
    
    // cout <<"-------------------------------------------"<<endl
    //     << "\n d1   = " << d1   << endl
    //     << "\n d2   = " << d1   << endl
    //     << "\n params.S   = " << d1   << endl
    //     << "\n params.K   = " << d1   << endl
    //     <<"-------------------------------------------"<<endl;

    return Ct;
}

//CRR Optimised Scaled Combinations - to force numbets to stay in scale
long double COSCer(long double N, long double j, long double quu, long double qdd) {
    if (j < 0 || j > N) return 0;

    // if (j > N - j) j = N - j;      // symmetry C(N,k) = C(N,N-k), I decrease iterations

    int um = 0, dm = 0;

    long double result = 1;
    for (long double i = 1; i <=j ; i++) {
        long double temu = 1.0, temd = 1.0;
        temu = quu;
        temd = qdd;
        um++;
        dm++;

        result *= ((N - (i - 1)) / i) * temd * temu;

    }
    result *= pow(qdd , (N - j) - j);
    dm += N - 2 * j;
    
    // dm+=int(1);
    // cout<<"\n - > COSCer for "<< N<<" = "<< j <<", "<< N-j <<", "<< um<<", "<< dm ; 

    return result;
}

// CDF of Normal Distribution using the error function
long double zScore(long double val){
    return 0.5 * erfc(-val * sqrt(1.0 / 2.0));
}


