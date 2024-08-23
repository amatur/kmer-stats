// using mpfr::mpreal for arbitrary precision arithmetic from https://github.com/advanpix/mpreal
// Compile with: g++ -std=c++11 -o thm12 thm12.cpp -lmpfr -lgmp

#include <iostream>
#include <random>
#include <string>
#include <set>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <map>
using namespace std;

#include "mpreal.h"
using mpfr::mpreal;
// Required precision of computations in decimal digits
const int digits = 50;  // Play with it to check different precisions

string type_map[] = {"random", "stress_random_lenk", "stress_ACGT", "stress_AC", "stress_A"};

//RESULTS
double estimator;
double error_term;
int type = 1;

namespace ValidationExamples{
    static const int NUM_EXAMPLES = 4;
    static const string strings[NUM_EXAMPLES] = {"AAAAAAAAA","ACGAACGAAA", "ACGTACGTAA","ACACACACAA"};
    static const string taus[NUM_EXAMPLES] = {"AAA", "AACG", "ACGT", "ACAA"};
    static const int ks[NUM_EXAMPLES] =  {3,4, 4,4};
    static const double rs[NUM_EXAMPLES] = {0.01, 0.01, 0.01, 0.05}; 
    static const double results_sds1000[NUM_EXAMPLES] = {0.770130, 0.193595, 0.275141, 0.413365};
};

class VarianceTest {
public:
    VarianceTest(std::string S, std::string tau, int k, double r) : S(S), tau(tau), k(k), r(r) {}
    VarianceTest(int EXAMPLE_ID) : EXAMPLE_ID(EXAMPLE_ID)  {}
    VarianceTest()  {}
    void doIt(int EXAMPLE_ID) {
        this->S = ValidationExamples::strings[EXAMPLE_ID];
        this->tau = ValidationExamples::taus[EXAMPLE_ID];
        this->k = ValidationExamples::ks[EXAMPLE_ID];
        this->r = ValidationExamples::rs[EXAMPLE_ID];
        this->EXAMPLE_ID = EXAMPLE_ID;
        doIt();
    }

    void doIt() {
        //int EXAMPLE_ID = 0;

        // Initialize the variables
        // this->S = ValidationExamples::strings[EXAMPLE_ID];
        // this->tau = ValidationExamples::taus[EXAMPLE_ID];
        // this->k = ValidationExamples::ks[EXAMPLE_ID];
        // this->r = ValidationExamples::rs[EXAMPLE_ID];
        int n = S.length();
        int L = n - k + 1;
        double p = r / (3 * (1 - r));
        double rr = pow(1 - r, k);

        double sum1 = 0;
        double sum2 = 0;
        double sum3 = 0;
        double temp = 0;

        for (int i = 0; i <= n - k; ++i) {
            int h1 = H(i);
            sum1 += (rr * pow(p, h1) - rr * rr * pow(p * p, h1));
            temp = 0;
            for (int j = i + 1; j <= i + k - 1; ++j) {
                if (j > n - k) continue;
                int h2 = H(j);
                sum2 += rr * rr * pow(p, h1) * pow(p, h2);
                temp += rr * rr * pow(p, h1) * pow(p, h2);

                int delta = k - (j - i);
                std::string over = overlap(delta);
                if (over.empty())
                    continue;
                int h3 = H(i, 2 * k - delta, over);
                sum3 += pow(1 - r, 2 * k - delta) * pow(p, h3);
            }
        }
        sum2 *= 2;
        sum3 *= 2;
        double sd = std::sqrt(sum1 - sum2 + sum3);
        std::printf("Computed: %.6f\n", std::sqrt(sum1 - sum2 + sum3));
        if(EXAMPLE_ID!=-1)
            cout<<"Expected: "<<ValidationExamples::results_sds1000[EXAMPLE_ID]<<endl;
        // std::printf("%.10f %.10f %.10f\n", sum1, sum2, sum3);
    }

private:
    std::string S;
    std::string tau;
    int k;
    double r;
    int EXAMPLE_ID = -1;

    std::string overlap(int delta) {
        bool o = true;
        for (int i = k - delta; i < k; ++i) {
            if (tau[i] != tau[i - k + delta]) {
                o = false;
                break;
            }
        }
        if (!o)
            return "";
        else
            return tau + tau.substr(delta, k - delta);
    }

    int H(int i, int l, const std::string& over) {
        int h = 0;
        for (int j = 0; j < l; ++j) {
            if (S[i + j] != over[j])
                ++h;
        }
        return h;
    }

    int H(int i) {
        int h = 0;
        for (int j = 0; j < k; ++j) {
            if (S[i + j] != tau[j])
                ++h;
        }
        return h;
    }
};


class SimParams{
    public:
        int k;
        unsigned long L;
        double r;

    SimParams(){
    }

    SimParams(int k, unsigned long L, double r){
        this->k = k;
        this->L = L;
        this->r = r;
    }
    double get_one_minus_q(){
        return pow(1.0-r, k);
        //(1 − r)^k
    }
    double get_q(){
        return 1 - pow(1.0-r, k);
        //1 − (1 − r)^k
    }
    double get_p(){
        return r/3.0/(1-r);
    }

    unsigned long get_n(){
        return L+k-1;
    }  
};

std::string generateRandomString(int length) {
    const std::string charset = "ACGT"; // Characters to choose from
    std::random_device rd;
    std::mt19937 gen(rd());

    //std::mt19937 gen(seed); to set seed
    std::uniform_int_distribution<> dis(0, charset.length() - 1);

    std::string randomString;
    randomString.reserve(length);

    for (int i = 0; i < length; ++i) {
        randomString += charset[dis(gen)];
    }

    return randomString;
}

std::string generateMutatedString(string s, double r) {
    string notT = "ACG";
    string notG = "ACT";
    string notC = "AGT";
    string notA = "CGT";

    std::random_device rd;
    std::mt19937 gen(rd());
    
    std::string mutatedString;
    mutatedString.reserve(s.length());

    for (int i = 0; i < s.length(); ++i) {
        std::uniform_real_distribution<> dis1(0.0, 1.0);
        // std::uniform_real_distribution<> dis2(0.0, 1.0);
        // std::uniform_real_distribution<> dis3(0.0, 1.0);
        double random_value = dis1(gen);
        
        
        //cout<<i<<" "<<random_value<<" "<<(1-r)<<endl;
        
        if (random_value < (1-r)) {
            mutatedString += s[i]; // Same with probability (1-r) //stay same
        } else {
            std::uniform_int_distribution<> dis(0, 2); // Range: {0, 1, 2}
            if(s[i] == 'A'){
                mutatedString += notA[dis(gen)];
            }else if(s[i] == 'C'){
                mutatedString += notC[dis(gen)];
            }else if(s[i] == 'G'){
                mutatedString += notG[dis(gen)];
            }else if(s[i] == 'T'){
                mutatedString += notT[dis(gen)];
            }
            
        }
    }
    return mutatedString;
}

/// Function to calculate the probability that a standard normal variable > c
double probability_greater_than(double c) {
    double cdf_c = 0.5 * (1 + std::erf(c / std::sqrt(2))); // Calculate the CDF of the standard normal distribution at c
    //cout<<"Z"<<1.0 - cdf_c<<endl; // The probability that a standard normal variable > c
    return 1.0 - cdf_c;
}

set<string> kspectrum(string s, int k){
    set<string> kmers;
    for(int i=0; i<s.length()-k+1; i++){
        string kmer = s.substr(i, k);
        kmers.insert(kmer);
    }
    return kmers;
}

int hammingDist(string str1, string str2) 
{ 
    int i = 0, count = 0; 
    while (str1[i] != '\0') { 
        if (str1[i] != str2[i]) 
            count++; 
        i++; 
    } 
    return count; 
} 

void tau_join(string kmer, vector<string>& retvec){
    int k = kmer.length();
    for (int delta = 1; delta<=k; delta++){
        bool isOverlap = true;
        for(int j = delta-1; j>=0; j--){
            if(kmer[j]!=kmer[k-delta+j]){
                isOverlap = false;
                break;
            }
        }
        if(isOverlap){
            if(delta==k){
                 retvec.push_back(kmer);
            }else{
                 retvec.push_back(kmer+kmer.substr(delta, k-delta));
            }
        }
    }
}

string generateStressTestString(int length, int k){
    std::string pattern = "ACGT";
    if(type == 1){
        pattern = generateRandomString(k); //say k=5, and randomstring is GTACTA: then GTACTAGTACTAGTACTA....
    }else if(type == 2){
        pattern = "ACGT";     //ACGTACGTACGT....
    }else if(type == 3){
        pattern = "AC";     //ACACACACACAC....
    }else if(type == 4){
        pattern = "A";   //AAAAAAAAAAAAA....
    }

    std::string s = "";
    std::string result;
    for (int i = 0; i < k; ++i) {
        result += pattern[i % pattern.length()];
    }

    for (int i = 0; i < floor(length/k); ++i) {
        s += result;
    }

    int patlen = s.length();
    for (int i = 1; i <= length-patlen; ++i) {
        s += "A";   // fill the remaining length with A
    }
    return s;
}

int intersect_size(string s, string t, int k){
    set<string> set1 = kspectrum(s, k);
    set<string> set2 = kspectrum(t, k);

    std::vector<string> intersect;

    std::set_intersection(set1.begin(), set1.end(),
                          set2.begin(), set2.end(),
                          std::back_inserter(intersect));
    return intersect.size();
}

string readString(string filename){
    string randomString;
    std::ifstream file(filename); 
    if (file.is_open()) {
        if (std::getline(file, randomString)) { 
            //std::cout << "Line read from file: " << line << std::endl;
            //L = randomString.length();
        } else {
            std::cerr << "File is empty." << std::endl;
        }
        file.close(); 
    } else {
        std::cerr << "Unable to open file." << std::endl;
    }
    return randomString;
}


/// @brief 
/// @param filename 
/// @param n 
/// @return string of length n from the file starting at a random position
string readSubstring(string filename, int n){
    std::ifstream file(filename);
    std::string line;

    if (file.is_open()) {
        file.seekg(0, std::ios::end); // Get the size of the file
        std::streampos fileSize = file.tellg();
        file.seekg(0, std::ios::beg);

        //std::srand(std::time(nullptr)); // Generate a random starting position within the file
        std::streampos startPos = (std::rand() % (int(fileSize) - n));
 
        // Move to the random starting position
        file.seekg(startPos);

        // Read 100 characters from the file
        char buffer[n+1]; // Buffer to hold 100 characters + null terminator
        file.read(buffer, n);
        buffer[n] = '\0'; // Null terminate the string

        // Assign the buffer content to the string
        line = buffer;

        //cout<< "sp" << startPos << " " <<std::rand() <<endl;
        //std::cout << "Random "<< L <<" characters from the file: " << line << std::endl;

        file.close(); // Close the file
    } else {
        std::cerr << "Unable to open file." << std::endl;
    }
    return line;
}


class ResultAggregator{
    public:
        vector<int> isizes;

        double isize_mean;
        double isize_sd;

        vector<double> test_vals_xtau;
        double test_vals_xtau_mean;
        double test_vals_xtau_sd;

        //vector<string> labels = {"r", "L","k", "num_reps", "estimate", "mean", "sd", "variance", "abs_error", "rel_error", "fixed_tau"};
        // = {"num_reps", "estimate", "mean", "sd", "variance", "abs_error", "rel_error"};

        string labels[9] = {"r", "L","k", "num_reps", "estimate", "mean", "sd", "variance", "fixed_tau"};
        vector<string> values;
        
        template<typename T>
        void put_values(T v){
            values.push_back(to_string(v));
        }

        double relative_error(double estimate, double mean){
            return abs((estimate-mean)/mean);
        }

        void diffReport(int num_replicates, double estimate, double mean, double sd, double r, int L, int k, string fixed_tau){
            put_values(r);
            put_values(L);
            put_values(k);
            put_values(num_replicates);
            put_values(estimate);
            put_values(mean);
            put_values(sd);
            put_values(sd*sd);
            // put_values(abs(estimate-mean));
            // put_values(relative_error(estimate, mean));
            values.push_back(fixed_tau);

            for(int i = 0; i< 9; i++){
                std::cout << std::fixed << std::setprecision(2);
                cout<<labels[i]<<" "<<values[i]<<" ";
            }
            cout<<endl;

            /*
            for(int i = 0; i< labels.size(); i++){
                cout<<values[i]<<" ";
            }
            cout<<endl;
            */
            
            //cout<<"num_reps "<< num_replicates<< " estimate "<<estimate<<" mean "<<mean<<" sd "<<sd<<" variance "<<sd*sd<<" abs_error "<<abs(estimate-mean)<<" rel_error "<<relative_error(estimate, mean)<<endl;
        }

        double calculateMean(const std::vector<int>& numbers) {
            double sum = 0.0;
            for (const auto& num : numbers) {
                sum += num;
            }
            return sum / numbers.size();
        }

        double calculateMean(const std::vector<double>& numbers) {
            double sum = 0.0;
            for (const auto& num : numbers) {
                sum += num;
            }
            return sum / numbers.size();
        }

        double calculateStandardDeviation(const std::vector<int>& numbers) {
            double mean = calculateMean(numbers);
            double squaredDiffSum = 0.0;
            for (const auto& num : numbers) {
                double diff = num - mean;
                squaredDiffSum += diff * diff;
            }
            double variance = squaredDiffSum / numbers.size();
            return std::sqrt(variance);
        }

        double calculateStandardDeviation(const std::vector<double>& numbers) {
            double mean = calculateMean(numbers);
            double squaredDiffSum = 0.0;
            for (const auto& num : numbers) {
                double diff = num - mean;
                squaredDiffSum += diff * diff;
            }
            double variance = squaredDiffSum / numbers.size();
            return std::sqrt(variance);
        }

        void calculateMeanSD(const std::vector<int>& numbers, double& mean, double& sd) {
            mean = calculateMean(numbers);
            double squaredDiffSum = 0.0;
            for (const auto& num : numbers) {
                double diff = num - mean;
                squaredDiffSum += diff * diff;
            }
            double variance = squaredDiffSum / numbers.size();
            sd = std::sqrt(variance);
        }

        void calculateMeanSD(const std::vector<double>& numbers, double& mean, double& sd) {
            mean = calculateMean(numbers);
            double squaredDiffSum = 0.0;
            for (const auto& num : numbers) {
                double diff = num - mean;
                squaredDiffSum += diff * diff;
            }
            double variance = squaredDiffSum / numbers.size();
            sd = std::sqrt(variance);
        }
};


class Thm11 {
    public:
        int k;
        double r;
        unsigned long L;
        int n;
        vector< map<int, string> > tau_joins_by_delta;
        int hd_s_tau[1000][1000];
        int unique = 0;

        Thm11() {}; // Provide a definition for the default constructor
        
        int get_hd_s_tau(int i, int j){
            if(i>=L || j>=unique){
                cout<<"Error: get_hd_s_tau out of bounds: i, j, L, unique: "<<i << " "<< j<<" "<<L<<" "<<unique<<endl;
                throw std::invalid_argument("Error: get_hd_s_tau out of bounds");
            }
            return hd_s_tau[i][j];
        }
        
        void set_hd_s_tau(int i, int j, int val){
            if(i>=L || j>=unique){
                cout<<"Error: get_hd_s_tau out of bounds: i, j, L, unique: "<<i << " "<< j<<" "<<L<<" "<<unique<<endl;
                throw std::invalid_argument("Error: set_hd_s_tau out of bounds");
            }
            hd_s_tau[i][j] = val;
        }

    void precompute_si_tau_hd_result_thm11(SimParams& sim, string& s){
        // // // PHASE1: INIT
        this->k = sim.k;
        this->r = sim.r;
        mpreal one_minus_r = 1 - r;        
        this->n = s.length();
        this->L = n - k + 1;
        mpreal p = r / (3 * (1 - r));
        mpreal rr = pow(one_minus_r, k);

        set<string> set1 = kspectrum(s, k);
        std::vector<string> s_vector(set1.begin(), set1.end());
        this->unique = s_vector.size();
        
        //populate tau_joins
        std::vector< vector<string> > tau_joins(s_vector.size());
        for(int j = 0 ; j < s_vector.size(); j++ ){
            string tau = s_vector[j]; //FiX 
            tau_join(tau, tau_joins[j]);
        }

        tau_joins_by_delta.resize(s_vector.size());
        for(int j = 0 ; j < s_vector.size(); j++ ){
            string tau = s_vector[j]; //FiX 
            for(string &s : tau_joins[j]){
                int delta = 2*k - s.length(); //2*k - s.length(); 
                // cout<<":::" << delta<<" "<<s<<endl;
                tau_joins_by_delta[j][delta] = s;
            }
        }
        
        // populate hd_s_tau
        for (int i = 0; i <= s.length() - k; ++i) {
            string s_i = (s.substr(i, k));
            for(int j = 0 ; j < s_vector.size(); j++ ){
                string tau = s_vector[j]; //FiX 
                set_hd_s_tau(i, j, hammingDist(s_i, tau));
            }
        }

        // // // PHASE2: Compute estimate of intersection
        mpreal sum_pr = 0;
        mpreal sum_e_max_tau = 0;

        for(int j = 0 ; j<s_vector.size(); j++){ // not on I, length of |sp(s)| -> so tau = s_vector[j]
            //string s_tau = s_vector[j];

            mpreal sum_top = 0;
            
            mpreal vcomp_sum1 = 0;
            mpreal vcomp_sum2 = 0;
            mpreal vcomp_sum3 = 0;
            mpreal vcomp_temp = 0;
            for (int i = 0; i <= n - k; ++i) {
                //PHASE1: do sum top
                mpreal pp = pow(p, get_hd_s_tau(i,j));
                mpreal pp2 = pow(pp,2);

                sum_top -= pp;

                // PHASE2: do sum bottom (sd)
                // double sum_bottom1 = 0;
                // double sum_bottom2 = 0;
                // double sum_bottom2_p1 = 0;
                // double sum_bottom2_p2 = 0;
                
                int h1 = get_hd_s_tau(i,j);
                vcomp_sum1 += (rr * pow(p, h1) - rr * rr * pow(p * p, h1));
                vcomp_temp = 0;
                for (int delta_o = i + 1; delta_o <= i + k - 1; ++delta_o) {
                    if (delta_o > n - k) continue; //fixing it : L-1 max value   n-k+1 -1 = n-k
                    int h2 = get_hd_s_tau(delta_o,j); //H(delta_o);
                    vcomp_sum2 += rr * rr * pow(p, h1) * pow(p, h2);
                    vcomp_temp += rr * rr * pow(p, h1) * pow(p, h2);

                    int delta = k - (delta_o - i);
                    std::string over="";
                    if(tau_joins_by_delta[j].count(delta)>0){
                        over = tau_joins_by_delta[j][delta]; // overlap(delta);
                    }
                    if (over=="")
                        continue;
                    int h3 = hammingDist(s.substr(i, 2 * k - delta), over); //H(i, 2 * k - delta, over); 

                    mpreal base_1minusr = 1 - r;
                    vcomp_sum3 += pow(base_1minusr, 2 * k - delta) * pow(p, h3);
                }
            }
            vcomp_sum2 *= 2;
            vcomp_sum3 *= 2;
            mpreal sd = sqrt(vcomp_sum1 - vcomp_sum2 + vcomp_sum3);
            //Finalized bottom for inner loop

            //cout<<":::"<<"sd: "<<sd<<" "<<s_vector[j]<<endl;
            

            sum_top *= rr;
            //Finalized top for inner loop

            sum_pr += probability_greater_than((sum_top / sd).toDouble());
            //sum_pr over all tau

            const double pi = M_PI; // Value of pi
            mpreal term1 = pow(2.0 / pi, 1.0 / 4.0);
            mpreal term2 = pow(2 * k, 2.0) * L / pow(sd, 3.0);
            mpreal term3 = sqrt(28.0) * pow(2 * k, 1.5) * sqrt(L) / (sqrt(pi) * pow(sd, 2.0));
            mpreal e_max_tau = term1 * sqrt(term2 + term3);
            //Finalized e_max_tau for inner loop

            sum_e_max_tau += e_max_tau;
        }
        estimator = sum_pr.toDouble();
        error_term = sum_e_max_tau.toDouble();
    }

    int sum_xi_tau(string& mutated_string, string tau){
        int sum = 0;
        // generate all kmer from mutated_string
        for (int i = 0; i <= mutated_string.length() - k; ++i) {
            std::string kmer = mutated_string.substr(i, k);
            if(kmer==tau){
                sum+=1;
            }
        }
        if(sum<0){
            cout<<"Error: overflow: sum_xi_tau "<<sum<<endl;
            exit(1);
        }
        return sum; //return CAPITAL-X^tau
    }

};


int main (int argc, char* argv[]){
    // Setup default precision for all subsequent computations
    // MPFR accepts precision in bits - so we do the conversion
    mpreal::set_default_prec(mpfr::digits2bits(digits));

    srand(time(nullptr));

    string filename="";
    int num_replicates = 1; //do this random string generation "num_replicates" times


    int L = 0;
    int k = 0;
    double r = 0;
    type = 0;
    num_replicates = 0;


    // /*** BEGIN EXAMPLE BYPASS COMMAND LINE ***/
    // Example: ./thm12 -l 100 -k 40 -r 0.01 -c 100
    // L=100;
    // k=40;
    // r=0.01;
    // num_replicates=100;
    // /*** END EXAMPLE BYPASS COMMAND LINE ***/
    
    vector<string> args(argv + 1, argv + argc);
    for (auto i = args.begin(); i != args.end(); ++i) {
        if (*i == "-h" || *i == "--help") {
            cout << "Syntax: -i [input-file] -l <length>[100] -r <mutation-rate>[0.01] -t <type-of-input-string>[1] -k <kmer-size>[33] -c <num-replicates>[100]" << endl;
            cout<< " Types: 0: random, 1: stress test, 2: ACGT, 3: AC, 4: A"<<endl;
            return 0;
        } else if (*i == "-r") {
            r = stod(*++i);
        } else if (*i == "-l") {
            L = stoull(*++i);
        } else if (*i == "-k") {
            k = stoi(*++i);
        } else if (*i == "-t") {
            type = stoi(*++i);
        }else if (*i == "-i") {
            filename = (*++i);
        } else if (*i == "-c") {
            num_replicates = stoi(*++i);
        }
    }

    // /*** BEGIN EXAMPLE BYPASS COMMAND LINE ***/
    // Example: ./thm12 -l 100 -k 40 -r 0.01 -c 100
    // L=100;
    // k=40;
    // r=0;
    // num_replicates=100;
    // /*** END EXAMPLE BYPASS COMMAND LINE ***/


    ResultAggregator res;
    std::string randomString;
    
    //Setup the FIXED string "s" and the parameters for the simulation (string in => randomString)
    SimParams sim(k, L, r); //k, L, r
    filename=""; //filename="/Users/amatur/code/downloads/t2t_chr21.sset";
    int length = sim.get_n();
    if (filename==""){
        if(type == 0){
            randomString = generateRandomString(length);
        }else{
            randomString = generateStressTestString(length, sim.k);
        }
    }else{
        randomString = readSubstring(filename, length);
    }

    // /*** BEGIN MINI TESTS ***/
    // int example_id = 2;
    // SimParams sim(ValidationExamples::ks[example_id], ValidationExamples::strings[example_id].length()-ValidationExamples::ks[example_id]+1, ValidationExamples::rs[example_id]); //k, L, r
    // randomString  = ValidationExamples::strings[example_id];
    // Thm11 thm11;
    // thm11.precompute_si_tau_hd_result_thm11(sim, randomString);
    // exit(1);
    // /*** END MINI TESTS ***/


    /*** BEGIN MAIN TEST ***/
    Thm11 thm11;
    thm11.precompute_si_tau_hd_result_thm11(sim, randomString);
    for(int i=0; i<num_replicates; i++){
        std::string mutatedString = generateMutatedString(randomString, sim.r);
        int Isize = intersect_size(randomString, mutatedString, sim.k);
        res.isizes.push_back(Isize);
    }
    res.calculateMeanSD(res.isizes, res.isize_mean, res.isize_sd);
    double r_prime = 1.0 - pow(2, log2( res.isize_mean)/sim.k-log2(sim.L)/sim.k);

    //cout.precision(digits);  
    bool NO_HEADER = false;
    if(NO_HEADER){
        cout<<num_replicates<< " "<< type_map[type] << " "<<sim.L<<" "<<sim.r<< " "<< sim.k <<  " " << estimator << " " << res.isize_mean << " " << res.isize_sd << " "<<sim.L*(1-sim.get_q())<< " "  << r_prime << " " << error_term << " "<< abs( res.isize_mean-estimator)<< " " << endl;
    }else{
        cout<<"num_replicates "<<num_replicates<< " input_s "<< type_map[type] << " L "<<sim.L<<" r "<<sim.r<< " k "<< sim.k <<  " estimate " << estimator << " |I|_mean " << res.isize_mean << " |I|_sd " << res.isize_sd << " L(1-q) "<<sim.L*(1-sim.get_q())<< " r_prime "  << r_prime << " theoretical_err " << error_term << " observed_err "<< abs( res.isize_mean-estimator)<< endl;
    }
}
