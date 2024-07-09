#include <iostream>
#include <random>
#include <string>
#include<set>
#include<fstream>

using namespace std;

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


set<string> kspectrum(string s, int k){
    set<string> kmers;
    for(unsigned long i=0; i<s.length()-k+1; i++){
        string kmer = s.substr(i, k);
        kmers.insert(kmer);
    }
    return kmers;
}

unsigned long intersect_size(string s, string t, int k){
    set<string> set1 = kspectrum(s, k);
    set<string> set2 = kspectrum(t, k);

    std::vector<string> intersect;

    std::set_intersection(set1.begin(), set1.end(),
                          set2.begin(), set2.end(),
                          std::back_inserter(intersect));
    return intersect.size();
}

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
    double get_q(){
        return 1 - pow(1.0-r, k);
        //1 − (1 − r)k
    }
    double get_p(){
        return r/3.0/(1-r);
    }

    unsigned long get_n(){
        return L+k-1;
    }

    
    
};

// k should be greater than alpha*log2(L)
double get_alpha(double r){
    return 1/(2 + log2(1/(1-r)));
}

class Thm3{

    public:
        int k;
        unsigned long L;
        double r;
        
    Thm3(SimParams p){
        this->k = p.k;
        this->L = p.L;
        this->r = p.r;
    }

    Thm3(int k, unsigned long L, double r){
        this->k = k;
        this->L = L;
        this->r = r;
    }


    // Overloading the assignment operator
    Thm3& operator=(const Thm3& other) {
        if (this != &other) {
            k = other.k;
            L = other.L;
            r = other.r;
        }
        return *this;
    }

    double error_term(){
        return (L - 1)*1.0/(pow(4.0*(1-r),k));
    }

    bool areEqual(double a, double b, double epsilon = 1e-9) {
        return fabs(a - b) < epsilon;
    }

    bool constraint_check(){
        //return areEqual(k/log2(L) - 1/(2 + log2(1/(1-r))) ,0);
        return (k/log2(L) > 1/(2 + log2(1/(1-r))));
    }
};

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

double calculateMean(const std::vector<int>& numbers) {
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


int main (int argc, char* argv[]){

    //do for different r
    srand(time(nullptr));

    unsigned long L = 10000;
    double alpha = -1;
    int k = -1;
    double r = 0.01;

    bool justGetMin = false;
    string filename="";

    vector<string> args(argv + 1, argv + argc);
    for (auto i = args.begin(); i != args.end(); ++i) {
            if (*i == "-h" || *i == "--help") {
                cout << "Syntax: -i <input-file> -a <alpha> -l <length> -r <mutation-rate> -p [getminalpha and exit] -k <kmer-size>" << endl;
                return 0;
            } else if (*i == "-a") {
                alpha = stod(*++i);
            } else if (*i == "-p") {
                justGetMin = true;
            } else if (*i == "-r") {
                r = stod(*++i);
            } else if (*i == "-l") {
                L = stoull(*++i);
            } else if (*i == "-k") {
                k = stoi(*++i);
            } else if (*i == "-i") {
                filename = (*++i);
            } 
    }

    if(justGetMin){
        double min_alpha = get_alpha(r);
        int min_k = log2(L)*min_alpha;
        cout<<"minAlpha "<< min_alpha<< " minK "<< min_k <<endl;
        return 0;
    }
    double min_alpha = get_alpha(r);
    if(alpha == -1){
        alpha = get_alpha(r);
    }
    if(k==-1){
        k = log2(L)*alpha;
    }
    
    Thm3 thm3(k, L, r);
    SimParams sim(k, L, r); //k, L, r
    int length = sim.get_n();

    std::string randomString;
    filename="";
    //filename="/Users/amatur/code/downloads/t2t_chr21.sset";
    if (filename!=""){
        randomString = readSubstring(filename, length);
    }
    double mean = 0;
    int num_replicates = 100; //do this random string generation 100 times
    vector<int> values(num_replicates);
    
    for(int i=0; i<num_replicates; i++){
        
        if (filename==""){
            randomString = generateRandomString(length);
        }
        //else{
        //     randomString = readSubstring(filename, length);
        // }
        
        std::string mutatedString = generateMutatedString(randomString, sim.r);
        int Isize = intersect_size(randomString, mutatedString, int(k));
        mean+=Isize;
        values[i] = Isize;
        //if(i<5)
        //cout<<int(k)<< "|I|: "<<intersect_size(randomString, mutatedString, int(k))<<endl;
    }
    mean /= 1.0*num_replicates;


    // std::cout << "Random string of length " << length << ": " << randomString << std::endl;
    // std::cout << "Mutated string of length " << length << ": " << mutatedString << std::endl;

    // cout<<"k"<< int(k)<< " |I|: "<<mean<<endl;
    // cout<<"L(1-q)"<<sim.L*(1-sim.get_q())<<endl;

    cout<<"r "<<r<<  " min_alpha " << min_alpha << " k "<< k << " alpha " << alpha<<" |I| "<<mean<<" L(1-q) "<<sim.L*(1-sim.get_q())<<" sd "<<calculateStandardDeviation(values)<<" err "<< thm3.error_term()<<" check "<<thm3.constraint_check()<<endl;

}