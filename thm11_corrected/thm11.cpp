#include <iostream>
#include <random>
#include <string>
#include<set>
#include<fstream>
#include<algorithm>
#include<iomanip>
#include<map>

using namespace std;
int hd_s_tau[1000][1000];

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
        double rr = std::pow(1 - r, k);

        double sum1 = 0;
        double sum2 = 0;
        double sum3 = 0;
        double temp = 0;

        for (int i = 0; i <= n - k; ++i) {
            int h1 = H(i);
            sum1 += (rr * std::pow(p, h1) - rr * rr * std::pow(p * p, h1));
            temp = 0;
            for (int j = i + 1; j <= i + k - 1; ++j) {
                if (j > n - k) continue;
                int h2 = H(j);
                sum2 += rr * rr * std::pow(p, h1) * std::pow(p, h2);
                temp += rr * rr * std::pow(p, h1) * std::pow(p, h2);

                int delta = k - (j - i);
                std::string over = overlap(delta);
                if (over.empty())
                    continue;
                int h3 = H(i, 2 * k - delta, over);
                sum3 += std::pow(1 - r, 2 * k - delta) * std::pow(p, h3);
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

class FastHammingDistance{
    public:
    string zodd="";
    string zeven="";
    uint64_t zodd_t, zeven_t;
    void hammingDistanceGlobalSet(int K){
        for (int i = 0; i < K+2; ++i) { //i<32 for K=30
            zodd += "01";
        }
        for (int i = 0; i < K+2; ++i) {
            zeven += "10";
        }
        zodd_t = std::stoull(zodd, nullptr, 2);
        zeven_t = std::stoull(zeven, nullptr, 2);
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


std::string encode_kmer(const std::string& kmer) {
    std::string encoded_kmer;
    for (char base : kmer) {
        switch (base) {
            case 'A':
                encoded_kmer += "00"; // 'A' is encoded as 00
                break;
            case 'C':
                encoded_kmer += "01"; // 'C' is encoded as 01
                break;
            case 'G':
                encoded_kmer += "10"; // 'G' is encoded as 10
                break;
            case 'T':
                encoded_kmer += "11"; // 'T' is encoded as 11
                break;
            default:
                throw std::invalid_argument("Invalid character in k-mer: " + std::string(1, base));
        }
    }
    return encoded_kmer;
}


vector<uint64_t> stringToBinaryVec(string kmer){
    int K = kmer.length();
    string bv_line = encode_kmer(kmer);
    vector<uint64_t> curr_bv(ceil(K/32));
    for (int i = 0; i< ceil(K/32); i++){
        if(2*K>i*64){
            curr_bv[i] = std::stoull(bv_line.substr(64*i,std::min(64, 2*K-i*64)), nullptr, 2);
        }
    }
    return curr_bv;
}

uint64_t stringToBinary(string kmer){
    string bv_line = encode_kmer(kmer);
    int K = kmer.length();
    uint64_t curr_bv_lo = std::stoull(bv_line.substr(0,std::min(64, 2*K)), nullptr, 2);
    return curr_bv_lo;
} 


std::string generateRandomString(int length, unsigned int seed=0) {
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

std::string generateMutatedString(string s, double r=0.1) {
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

// Function to calculate the probability that a standard normal variable > c
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
    
    for (int delta = 1; delta<=k; delta++){ // <= k 
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
            //len=|K|+|K|-delta
            //
        }
       
        // kmer[0] == kmer[k-1]

        // kmer[0] = kmer[k-2]
        // kmer[1] = kmer[k-1]

        // kmer[0] = kmer[k-3]
        // kmer[1] = kmer[k-2]
        // kmer[2] = kmer[k-1]
    }
    
}

double pow(double base, int exp){
    double val = pow(base, (double)exp);
    //cout<< base<< "^" << exp<<"=" <<val<<endl; 

    if(abs(base)>0 && abs(base)<1 && exp>0 && val > base){
        cout<<"error pow "<<val<<endl;
        throw ("error");
    }

    if(val != 0.0 && val < std::numeric_limits<double>::min()){
        cout<<"underflow "<<val<<endl;
        throw std::underflow_error("underflow");
    }

    return val;
}
void value_checker(double val){
    if(val != 0.0 && val < std::numeric_limits<double>::min()){
        cout<<"underflow "<<val<<endl;
        throw std::underflow_error("underflow");
    }
    if(std::isinf(val)){
        cout<<"overflow "<<val<<endl;
        throw std::overflow_error("overflow");
    }
}



string generateStressTestString(int length, int k){
    std::string pattern = "ACGT";
    if(type == 1){
        pattern = generateRandomString(k);
    }else if(type == 2){
        pattern = "ACGT";
    }else if(type == 3){
        pattern = "AC";
    }else if(type == 4){
        pattern = "A";
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
        s += "A";
    }
    //cout<<length<<" " << s.length()<<" "<<"stress test str"<<s<<endl;
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
string readSubstring(string filename, int L){
    std::ifstream file(filename);
    std::string line;

    if (file.is_open()) {
        file.seekg(0, std::ios::end); // Get the size of the file
        std::streampos fileSize = file.tellg();
        file.seekg(0, std::ios::beg);

        //std::srand(std::time(nullptr)); // Generate a random starting position within the file
        std::streampos startPos = (std::rand() % (int(fileSize) - L));
 
        // Move to the random starting position
        file.seekg(startPos);

        // Read 100 characters from the file
        char buffer[L+1]; // Buffer to hold 100 characters + null terminator
        file.read(buffer, L);
        buffer[L] = '\0'; // Null terminate the string

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
        vector<string> labels = {"r", "L","k", "num_reps", "estimate", "mean", "sd", "variance", "fixed_tau"};

        vector<string> values;

    
        template<typename T>
        void put_values(T v){
            values.push_back(to_string(v));
        }


        // = {"num_reps", "estimate", "mean", "sd", "variance", "abs_error", "rel_error"};


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

        for(int i = 0; i< labels.size(); i++){
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
        Thm11() {}; // Provide a definition for the default constructor

        int k = 8;
        double r = 0.01;
        unsigned long L = 17;
        int n = L+k-1;
        int fixed_tau_index_j = 0;
        
        //set by call to init
        string fixed_string_s; 
        SimParams sim;
        
        vector< map<int, string> > tau_joins_by_delta;

        int hd_s_tau[1000][1000];

        string fixed_tau;
        int unique = 0;

        //saved 1 : string: TACGTACGAA, tau ACGA
        
        int get_hd_s_tau(int i, int j){
            if(i>=n || j>=unique){
                cout<<"Error: get_hs_s_tau out of bounds"<<endl;
                exit(1);
            }
            return hd_s_tau[i][j];
        }
        
    void precompute_si_tau_hd_result_thm11(SimParams& sim, string& s){
        int k = sim.k;
        double r = sim.r;
        // double p = sim.get_p();
        // int L = s.length() - k + 1;
        
        int n = s.length();
        int L = n - k + 1;
        double p = r / (3 * (1 - r));
        double rr = std::pow(1 - r, k);

        set<string> set1 = kspectrum(s, k);
        std::vector<string> s_vector(set1.begin(), set1.end());
        unique = s_vector.size();
        
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
                cout<<delta<<" "<<s<<endl;
                tau_joins_by_delta[j][delta] = s;
            }
        }
        
        //populate hd_s_tau
        for (int i = 0; i <= s.length() - k; ++i) {
            string s_i = (s.substr(i, k));
            for(int j = 0 ; j < s_vector.size(); j++ ){
                string tau = s_vector[j]; //FiX 
                hd_s_tau[i][j] = hammingDist(s_i, tau);
            }
        }


        double sum_pr = 0;
        double sum_e_max_tau = 0;

        for(int j = 0 ; j<s_vector.size(); j++){ // not on I, length of |sp(s)| -> so tau = s_vector[j]
            //string s_tau = s_vector[j];

            double sum_top = 0;
            
            double vcomp_sum1 = 0;
            double vcomp_sum2 = 0;
            double vcomp_sum3 = 0;
            double vcomp_temp = 0;
            for (int i = 0; i <= n - k; ++i) {
                //PHASE1: do sum top
                double pp = pow(p, get_hd_s_tau(i,j));
                double pp2 = pow(pp,2);

                sum_top -= pp;

                // PHASE2: do sum bottom (sd)
                // double sum_bottom1 = 0;
                // double sum_bottom2 = 0;
                // double sum_bottom2_p1 = 0;
                // double sum_bottom2_p2 = 0;
                
                int h1 = get_hd_s_tau(i,j);
                vcomp_sum1 += (rr * std::pow(p, h1) - rr * rr * std::pow(p * p, h1));
                vcomp_temp = 0;
                for (int delta_o = i + 1; delta_o <= i + k - 1; ++delta_o) {
                    if (delta_o > n - k) continue;
                    int h2 = get_hd_s_tau(delta_o,j); //H(delta_o);
                    vcomp_sum2 += rr * rr * std::pow(p, h1) * std::pow(p, h2);
                    vcomp_temp += rr * rr * std::pow(p, h1) * std::pow(p, h2);

                    int delta = k - (delta_o - i);
                    std::string over="";
                    if(tau_joins_by_delta[j].count(delta)>0){
                        over = tau_joins_by_delta[j][delta]; // overlap(delta);
                    }
                    if (over=="")
                        continue;
                    int h3 = hammingDist(s.substr(i, 2 * k - delta), over); //H(i, 2 * k - delta, over); 
                    vcomp_sum3 += std::pow(1 - r, 2 * k - delta) * std::pow(p, h3);
                }
            }
            vcomp_sum2 *= 2;
            vcomp_sum3 *= 2;
            double sd = std::sqrt(vcomp_sum1 - vcomp_sum2 + vcomp_sum3);
            cout<<"sd: "<<sd<<endl;
            //Finalized bottom

            sum_top *= rr;
            //Finalized top

            //sum_pr += probability_greater_than(sum_top / sqrt(sum_bottom1 + sum_bottom2));
            sum_pr += probability_greater_than(sum_top / sd);
            //sum_pr over all tau


            const double pi = M_PI; // Value of pi
            double term1 = pow(2.0 / pi, 1.0 / 4.0);
            double term2 = pow(2 * k, 2.0) * L / pow(sd, 3.0);
            double term3 = sqrt(28.0) * pow(2 * k, 1.5) * sqrt(L) / (sqrt(pi) * pow(sd, 2.0));
            double e_max_tau = term1 * sqrt(term2 + term3);
            //Finalized e_max_tau

            //cout<<term1<<" "<<term2<<" "<<term3<<" "<<e_max_tau<<endl;
            sum_e_max_tau += e_max_tau;
            //sum_e_max_tau over all tau

            /*
            for(int i = 0; i< L; i++){ 
                double pp=pow(p, hd_s_tau[i][j]);
                double pp2 = pow(pp,2);

                sum_top -= pp;
                sum_bottom1 += pp/sim.get_one_minus_q() - pp2;
                //sum_bottom1 += pp/();


                for(int delta = 1; delta<=k; delta++){
                    // if(i+delta>=L){
                    //     break;
                    // }
                    sum_bottom2_p1+= -2*pow(p, hd_s_tau[i][j] + hd_s_tau[i+k-delta][j] );
                    //cout<<sum_bottom2_p1<<" sp "<< " delta " << delta << " k "<<k<< " "<< hd_s_tau[i][j] << " "<< hd_s_tau[i+delta][j]  << " " << i+delta<<endl;
                }
                
                for (int t = 0; t< tau_joins[j].size(); t++){
                    string joinedString = tau_joins[j][t];
                    int delta = 2*k - joinedString.length();//len = 2k -delta

                    //sum_bottom2+= -2*pow(p, hd_s_tau[i][j] + hd_s_tau[i+delta][j] );

                    string s_i_len_same_as_joinedString = s.substr(i, joinedString.length() );
                    if(s_i_len_same_as_joinedString.length()==joinedString.length() && delta!=k){
                        int d = hammingDist(s_i_len_same_as_joinedString, joinedString);
                        sum_bottom2_p2 += 2*pow(p,d)/(pow(1-r, delta));
                    }

                    
                }
            }
            sum_bottom2= sum_bottom2_p1 + sum_bottom2_p2;

            cout<<"for tau "<<j<<" "<<s_vector[j]<<endl;
            cout<<"sum top = "<<sum_top<<endl;
            cout<<"sum b1 = "<<sum_bottom1<<endl;
            cout<<"sum b2 p1 = "<<sum_bottom2_p1<<endl;
            cout<<"sum b2 p2 = "<<sum_bottom2_p2<<endl;
            

            cout<<"sum b2 ="<<sum_bottom2<<endl;

            cout<<"sum_top / sqrt(sum_bottom1 + sum_bottom2 = "<<sum_top / sqrt(sum_bottom1 + sum_bottom2)<<endl;



            if (std::isinf(sum_bottom2_p1) && sum_bottom2_p1 < 0) {
                cout<<"Error "<<j<<" "<< s_vector[j]<<endl;
                exit(1);
            }
            double sigma_tau = sqrt(sum_bottom1 + sum_bottom2); //variance of tau

            const double pi = M_PI; // Value of pi
            double term1 = pow(2.0 / pi, 1.0 / 4.0);
            double term2 = pow(2 * k, 2.0) * L / pow(sigma_tau, 3.0);
            double term3 = sqrt(28.0) * pow(2 * k, 1.5) * sqrt(L) / (sqrt(pi) * pow(sigma_tau, 2.0));
            double e_max_tau = term1 * sqrt(term2 + term3);
            //cout<<term1<<" "<<term2<<" "<<term3<<" "<<e_max_tau<<endl;
            sum_e_max_tau += e_max_tau;
            sum_pr += probability_greater_than(sum_top / sqrt(sum_bottom1 + sum_bottom2));
            */
            
            
        }
        estimator = sum_pr;
        error_term = sum_e_max_tau;

        //cout<<"Sum Pr"<<sum_pr<<"I "<<I<<endl;
    }

    double compute_sd_estimate(){
        //PART1:: INIT
        // fixed_string_s = readSubstring("/Users/amatur/code/downloads/t2t_chr21.sset", L+k-1);
        //fixed_string_s = "AGCTTAAAGTAATTATCTAGGTGTCTGTATTTG";
        //fixed_string_s = "AGCTTAAAGTAATTATCTAGGTGTCTGTATTTGCCT";
        vector<string> s_kspectrum_vector;
    // L = fixed_string_s.length() - k + 1;
        //readSubstring("/Users/amatur/code/downloads/t2t_chr21.sset", L+k-1);
        type = 3;
        //fixed_string_s = generateStressTestString( L+k-1, k);
        fixed_string_s = generateRandomString( L+k-1);
        //generateStressTestString( L+k-1, k);
        
        //fixed_string_s = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG";
        //  fixed_string_s = "ACGAACGAAA";
        //  L = fixed_string_s.length() - k + 1;

        fixed_string_s = "AAAAAAAAA";


        cout<<fixed_string_s<<endl;
        sim.k = k;
        sim.L = L;
        sim.r = r;
        set<string> set1 = kspectrum(fixed_string_s, k);
        std::vector<string> s_vector(set1.begin(), set1.end());
        s_kspectrum_vector = s_vector;

        
        fixed_tau = s_kspectrum_vector[fixed_tau_index_j];
        if(fixed_tau_index_j > s_kspectrum_vector.size()){
            cout<<"Error: fixed_tau_index_j out of bounds"<<" " <<fixed_tau_index_j<<endl;
            exit(1);
        }

        vector< vector<string> > tau_joins;
        tau_joins.resize(s_kspectrum_vector.size());
        for(int j = 0 ; j < s_kspectrum_vector.size(); j++ ){
            string tau = s_kspectrum_vector[j]; //FiX 
            tau_join(tau, tau_joins[j]);
        }
        for(int j = 0 ; j < s_kspectrum_vector.size(); j++ ){
            string tau = s_kspectrum_vector[j]; //FiX 
            for(string &s : tau_joins[j]){
                int delta = 2*k - s.length(); //2*k - s.length(); 
                tau_joins_by_delta[j][delta] = s;
            }
        }

        for (int i = 0; i < fixed_string_s.size(); ++i) {
            string s_i = (fixed_string_s.substr(i, k));
            //for(int j = 0 ; j < s_kspectrum_vector.size(); j++ ){
            for(int j = fixed_tau_index_j ; j <= fixed_tau_index_j; j++ ){
                string tau = s_kspectrum_vector[j]; //FiX 
                hd_s_tau[i][j] = hammingDist(s_i, tau);
                //tau_join(tau, tau_joins[j]);
            }
        }

        
        //PART2:: COMPUTE
        double p = sim.get_p();
        double sigma_tau = 777777;
        //fix a tau: j = 0
        for(int j = fixed_tau_index_j ; j<=fixed_tau_index_j; j++){ // not on I, length of |sp(s)|
            double sum_top = 0;
            double sum_bottom1 = 0;
            double sum_bottom2 = 0;
            double sum_bottom2_p1 = 0;
            double sum_bottom2_p2 = 0;

            set<int> delta_taus;
             for (int t = 0; t< tau_joins[j].size(); t++){
                  delta_taus.insert(2*k - tau_joins[j][t].length());
             }

            int coutcounter=0;
            for(int i = 0; i< L; i++){ 
                double pp = pow(p, hd_s_tau[i][j]);
                double pp2 = pow(pp,2);

                sum_top -= pp;
                sum_bottom1 += pp/sim.get_one_minus_q() - pp2;
                //sum_bottom1 += pp/();

                int t =  tau_joins[j].size()-1;
                string joinedString = tau_joins[j][t];
                int delta = 2*k - joinedString.length();//len = 2k -delta

                for(int o = 1; o <= k-1; o++){
                    if(L-i-1 < o){
                        cout<<"SHOUDL NOT HAPPEN"<<endl;
                        break;
                    }
                    sum_bottom2_p1+= -2*pow(p, hd_s_tau[i][j] + hd_s_tau[i+k-o][j] );
                    //cout<<"CC" << p<<" "<< fixed_string_s.substr(i, k) << " "<< o << " hd "<< hd_s_tau[i][j] << " " << fixed_string_s.substr(i+k-o, k) <<" hd "<< hd_s_tau[i+k-o][j] <<" " << -2*pow(p, hd_s_tau[i][j] + hd_s_tau[i+k-o][j] )<<endl;
                    
                    //cout<<sum_bottom2_p1<<" sp "<< " delta " << delta << " k "<<k<< " "<< hd_s_tau[i][j] << " "<< hd_s_tau[i+delta][j]  << " " << i+delta<<endl;
                    //if(i>L+k-1 = o upper range = L - i )
                    
                    /*
                    if(k-o == delta && t>=0){
                        string s_i_len_same_as_joinedString = fixed_string_s.substr(i, joinedString.length() );
                        if(s_i_len_same_as_joinedString.length()==joinedString.length() && delta!=k){
                            coutcounter++;
                            int d = hammingDist(s_i_len_same_as_joinedString, joinedString);
                            //cout<<d<< " " << i<<" "<< s_i_len_same_as_joinedString << " "<< delta << " hd "<< hd_s_tau[i][j] << " " << joinedString<<" hd "<< hd_s_tau[i+k-delta][j] <<" " << -2*pow(p, hd_s_tau[i][j] + hd_s_tau[i+k-delta][j] )<<endl;

                            sum_bottom2_p2 += 2*pow(p,d)*(pow(1-r, k+o));
                            //cout<<(pow(1-r, delta))<<endl;
                        }
                                                //problem in p2
                        //go to next delta
                        t--;
                        if(t>=0){
                            joinedString = tau_joins[j][t];
                            delta = 2*k - joinedString.length();//len = 2k -delta
                        }
                    }
                    */
                    //sum_bottom2+= -2*pow(p, hd_s_tau[i][j] + hd_s_tau[i+delta][j] );
                    //sum_bottom2_p2 += 2*pow(p,d)/(pow(1-r, delta));
                }
                
                ///*
                for (int t = 0; t< tau_joins[j].size(); t++){
                    string joinedString = tau_joins[j][t];
                    int delta = 2*k - joinedString.length();//len = 2k -delta
                    if (delta==k)
                    {
                        continue;
                    }
                    

                    //sum_bottom2+= -2*pow(p, hd_s_tau[i][j] + hd_s_tau[i+delta][j] );
                    string s_i_len_same_as_joinedString  ="";
                    if(joinedString.length() - i + 1 > 0 && i<fixed_string_s.length()){
                        s_i_len_same_as_joinedString = fixed_string_s.substr(i, joinedString.length() );
                    }

                    if(s_i_len_same_as_joinedString.length()==joinedString.length()){
                    //if(delta!=k){
                        coutcounter++;
                        int d = hammingDist(s_i_len_same_as_joinedString, joinedString);
                        //cout<<d<< " " << i<<" "<< s_i_len_same_as_joinedString << " "<< delta << " hd "<< hd_s_tau[i][j] << " " << joinedString<<" hd "<< hd_s_tau[i+k-delta][j] <<" " << -2*pow(p, hd_s_tau[i][j] + hd_s_tau[i+k-delta][j] )<<endl;

                        sum_bottom2_p2 += 2*pow(p,d)/(pow(1-r, delta));
                        //cout<<(pow(1-r, delta))<<endl;
                    }
                    //problem in p2
                    
                }
                //*/
            }
            cout<<"Count "<<coutcounter<<endl;
            sum_bottom2= sum_bottom2_p1 + sum_bottom2_p2;
            cout<<"sum top "<<sum_top<<endl;
            cout<<"sum b1 (p0) "<<sum_bottom1 / pow((1-r), 2*k) <<endl;
            cout<<"sum b2 p1 "<<sum_bottom2_p1 / pow((1-r), 2*k) <<endl;
            cout<<"sum b2 p2 "<<sum_bottom2_p2/ pow((1-r), 2*k) <<endl;
            double a = sum_bottom1 / pow((1-r), 2*k) ;//0.201732;
            double b = sum_bottom2_p1  / pow((1-r), 2*k);
            //-20.712563;
            //sum_bottom2_p1 / pow((1-r), 2*k) ;
            //-20.712563;
            double c = sum_bottom2_p2/ pow((1-r), 2*k);
            //21.037053;
            //

            cout<<"sum b2 "<<sum_bottom2<<endl;

            cout<<"sum_top / sqrt(sum_bottom1 + sum_bottom2) "<<sum_top / sqrt(sum_bottom1 + sum_bottom2)<<endl;

            if (std::isinf(sum_bottom2_p1) && sum_bottom2_p1 < 0) {
                cout<<"Error "<<j<<" "<< s_kspectrum_vector[j]<<endl;
                exit(1);
            }
            //sigma_tau = pow((1-r), k)  * sqrt(sum_bottom1 + sum_bottom2); //variance of tau
            sigma_tau = pow((1-r), k)  * sqrt(sum_bottom1 + sum_bottom2); //variance of tau
            cout<<"ABC"<<sqrt(a+b+c)<<endl;
            //sigma_tau = sqrt(a+b+c);
            
        }
        return sigma_tau;
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



///*
int main (int argc, char* argv[]){

    // VarianceTest varianceTest;
    // varianceTest.doIt(0);
    // exit(1);


    //do for different r
    srand(time(nullptr));

    int L = -1;
    int k = -1;
    double r = -1;

    string filename="";


    int num_replicates = 1; //do this random string generation 100 times
    vector<string> args(argv + 1, argv + argc);
    for (auto i = args.begin(); i != args.end(); ++i) {
            if (*i == "-h" || *i == "--help") {
                cout << "Syntax: -i [input-file] -l <length> -r <mutation-rate> -t <type-of-input-string> -k <kmer-size> -c <num-replicates>" << endl;
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

    //hammingDistanceGlobalSet(K);
    

    // L=100;
    // k=40;
    // r=0;
    // SimParams sim(k, L, r); //k, L, r
    num_replicates=100;
    
    
    ResultAggregator res;
    std::string randomString;
    
    //do a single string
    /*
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
    */

    int example_id = 2;
    SimParams sim(ValidationExamples::ks[example_id], ValidationExamples::strings[example_id].length()-ValidationExamples::ks[example_id]+1, ValidationExamples::rs[example_id]); //k, L, r
    randomString  = ValidationExamples::strings[example_id];
    Thm11 thm11;
    thm11.precompute_si_tau_hd_result_thm11(sim, randomString);
    //exit(1);

    for(int i=0; i<num_replicates; i++){
        std::string mutatedString = generateMutatedString(randomString, sim.r);
        int Isize = intersect_size(randomString, mutatedString, sim.k);
        res.isizes.push_back(Isize);
    }

    res.calculateMeanSD(res.isizes, res.isize_mean, res.isize_sd);
    
    double r_prime = 1.0 - pow(2, log2( res.isize_mean)/sim.k-log2(sim.L)/sim.k);
    cout<<"num_replicates "<<num_replicates<< " input_s "<< type_map[type] << " L "<<sim.L<<" r "<<sim.r<< " k "<< sim.k <<  " estimate " << estimator << " |I|_mean " << res.isize_mean << " |I|_sd " << res.isize_sd << " L(1-q) "<<sim.L*(1-sim.get_q())<< " r_prime "  << r_prime << " theoretical_err " << error_term << " observed_err "<< abs( res.isize_mean-estimator)<< endl;
    //cout<<endl;

}
//*/