#include <iostream>
#include <random>
#include <string>
#include<set>




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


int main() {
    
    
    double alpha = 10;
    unsigned long L = 10000;
    double k = log2(1000)*alpha;
    double r = 0.01;


    Thm3 thm3(k, L, r);
    cout<<thm3.k<<" "<<thm3.L<<" " << thm3.error_term()<<" "<<thm3.constraint_check()<<endl;

    
    k = log2(L)*alpha;
    thm3 = Thm3(k, L, 0.01);
    cout<<thm3.k<<" "<<thm3.L<<" "<<thm3.error_term()<<" "<<thm3.constraint_check()<<endl;


    alpha = get_alpha(r);
    cout<<"Min Alpha "<< alpha<<endl;//r = 0.01 alpha = 0.5
    alpha = 70; //alpha 50
    cout<<"Alpha "<< alpha<<endl;
    k = log2(L)*alpha;
    thm3 = Thm3(k, L, r);
    cout<<thm3.k<<" "<<thm3.L<<" "<<thm3.error_term()<<" "<<thm3.constraint_check()<<endl;


    SimParams sim(k, L, r); //k, L, r
    int length = sim.get_n();

    

    std::string randomString = generateRandomString(length);
    std::string mutatedString = generateMutatedString(randomString, sim.r);
    //do this random string generation 1000 times
    double mean = 0;

    int num_replicates = 100;
    for(int i=0; i<num_replicates; i++){
        randomString = generateRandomString(length);
        mutatedString = generateMutatedString(randomString, sim.r);
        mean+=intersect_size(randomString, mutatedString, int(k));

        if(i<5)
        cout<<int(k)<< "|I|: "<<intersect_size(randomString, mutatedString, int(k))<<endl;
        //cout<<"L(1-q)"<<sim.L*(1-sim.get_q())<<endl;
    }
    mean/=1.0*num_replicates;


    // std::cout << "Random string of length " << length << ": " << randomString << std::endl;
    // std::cout << "Mutated string of length " << length << ": " << mutatedString << std::endl;

    cout<<"k"<< int(k)<< " |I|: "<<mean<<endl;
    cout<<"L(1-q)"<<sim.L*(1-sim.get_q())<<endl;



    // TODO
    // do a plot alpha, k, |I|, L(1-q)
    // run it on ecoli 10000 bp
}