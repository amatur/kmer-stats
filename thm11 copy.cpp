#include <iostream>
#include <random>
#include <string>
#include<set>
#include<fstream>

using namespace std;
//std::string randomString;
int hd_si_tau[1000][1000];


string zodd="";
string zeven="";
uint64_t zodd_t, zeven_t;
string type_map[] = {"random", "stress_random_lenk", "stress_ACGT", "stress_AC", "stress_A"};


//RESULTS
double estimator;
int S_obs;
double error_term;
int type = 1;


class SimParams{
    public:
    int k;
    unsigned long L;
    double r;

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
// function to calculate Hamming distance 
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
            retvec.push_back(kmer+kmer.substr(delta, k-delta));
        }
       
        // kmer[0] == kmer[k-1]

        // kmer[0] = kmer[k-2]
        // kmer[1] = kmer[k-1]

        // kmer[0] = kmer[k-3]
        // kmer[1] = kmer[k-2]
        // kmer[2] = kmer[k-1]
    }
    
}

void precompute_si_tau_hd_result_thm11(SimParams& sim, string& s){
    int k = sim.k;
    set<string> set1 = kspectrum(s, k);
    std::vector<string> s_vector(set1.begin(), set1.end());

    std::vector< vector<string> > tau_joins(s_vector.size());
    for (int i = 0; i <= s.length() - k; ++i) {
        string s_i = (s.substr(i, k));
        for(int j = 0 ; j < s_vector.size(); j++ ){
            string tau = s_vector[j]; //FiX 
            hd_si_tau[i][j] = hammingDist(s_i, tau);
            tau_join(tau, tau_joins[j]);
        }
    }

    double p = sim.get_p();
    int L = s.length() - k + 1;
    double r = sim.r;

    double sum_pr = 0;
    double sum_e_max_tau = 0;
    for(int j = 0 ; j<s_vector.size(); j++){ // not on I, length of |sp(s)|

        double sum_top = 0;
        double sum_bottom1 = 0;
        double sum_bottom2 = 0;

        
        for(int i = 0; i< L; i++){ 
            double pp=pow(p, hd_si_tau[i][j]);
            double pp2 = pow(pp,2);

            sum_top -= pp;
            sum_bottom1 += pp/sim.get_one_minus_q() - pp2;
            //sum_bottom1 += pp/();


            for(int delta = 1; delta<=k; delta++){
                sum_bottom2+= -2*pow(p, hd_si_tau[i][j] + hd_si_tau[i+delta][j] );
            }

            for (int t = 0; t< tau_joins[j].size(); t++){
                string joinedString = tau_joins[j][t];
                int delta = 2*k - joinedString.length();//len = 2k -delta

                //sum_bottom2+= -2*pow(p, hd_si_tau[i][j] + hd_si_tau[i+delta][j] );

                string s_i_len_same_as_joinedString = s.substr(i, joinedString.length() );
                int d = hammingDist(s_i_len_same_as_joinedString, joinedString);
                sum_bottom2 += 2*pow(p,d)/(pow(1-r, delta));
            }
        }
        // cout<<"sum top "<<sum_top<<endl;
        // cout<<"sum b1 "<<sum_bottom1<<endl;
        // cout<<"sum b2 "<<sum_bottom2<<endl;
        // cout<<"sum_top / sqrt(sum_bottom1 + sum_bottom2 "<<sum_top / sqrt(sum_bottom1 + sum_bottom2)<<endl;

        double sigma_tau = sqrt(sum_bottom1 + sum_bottom2);
        
        
        const double pi = M_PI; // Value of pi
        double term1 = pow(2.0 / pi, 1.0 / 4.0);
        double term2 = pow(2 * k, 2.0) * L / pow(sigma_tau, 3.0);
        double term3 = sqrt(28.0) * pow(2 * k, 1.5) * sqrt(L) / (sqrt(pi) * pow(sigma_tau, 2.0));
        double e_max_tau = term1 * sqrt(term2 + term3);
        //cout<<term1<<" "<<term2<<" "<<term3<<" "<<e_max_tau<<endl;
        sum_e_max_tau += e_max_tau;
        sum_pr += probability_greater_than(sum_top / sqrt(sum_bottom1 + sum_bottom2));
    }
    estimator = sum_pr;
    error_term = sum_e_max_tau;

    //cout<<"Sum Pr"<<sum_pr<<"I "<<I<<endl;
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

// double calculateMean(const std::vector<int>& numbers) {
//     double sum = 0.0;
//     for (const auto& num : numbers) {
//         sum += num;
//     }
//     return sum / numbers.size();
// }

// double calculateStandardDeviation(const std::vector<int>& numbers) {
//     double mean = calculateMean(numbers);
//     double squaredDiffSum = 0.0;
//     for (const auto& num : numbers) {
//         double diff = num - mean;
//         squaredDiffSum += diff * diff;
//     }
//     double variance = squaredDiffSum / numbers.size();
//     return std::sqrt(variance);
// }

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

class ResultAggregator{
    public:
        vector<int> isizes;

    double isize_mean;
    double isize_sd;

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


int main (int argc, char* argv[]){
    //do for different r
    srand(time(nullptr));

    int L = -1;
    int k = -1;
    double r = -1;

    string filename="";

    vector<string> args(argv + 1, argv + argc);
    for (auto i = args.begin(); i != args.end(); ++i) {
            if (*i == "-h" || *i == "--help") {
                cout << "Syntax: -i [input-file] -l <length> -r <mutation-rate> -t <type-of-input-string> -k <kmer-size>" << endl;
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
            } 
    }

    //hammingDistanceGlobalSet(K);
    
    SimParams sim(5, 100, 0.01); //k, L, r
    //SimParams sim(k, L, r); //k, L, r

    int num_replicates = 100; //do this random string generation 100 times
    ResultAggregator res;

    //do a single string
    std::string randomString;
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
    
    for(int i=0; i<num_replicates; i++){
        std::string mutatedString = generateMutatedString(randomString, sim.r);
        precompute_si_tau_hd_result_thm11(sim, randomString);
        int Isize = intersect_size(randomString, mutatedString, sim.k);
        res.isizes.push_back(Isize);
    }

    res.calculateMeanSD(res.isizes, res.isize_mean, res.isize_sd);
    
    double r_prime = 1.0 - pow(2, log2(S_obs)/sim.k-log2(sim.L)/sim.k);
    cout<<"input_s "<< type_map[type] << " L "<<sim.L<<" r "<<sim.r<< " k "<< sim.k <<  " estimate " << estimator << " |I|_mean " << res.isize_mean << " |I|_sd " << res.isize_sd << " L(1-q) "<<sim.L*(1-sim.get_q())<< " r_prime "  << r_prime << " theoretical_err " << error_term << " observed_err "<< abs( res.isize_mean-estimator)<< endl;
    //cout<<endl;

}