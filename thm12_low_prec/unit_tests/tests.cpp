
#include <iostream>
#include <string>
#include <cmath>
#include<iomanip>
#include "thm11.cpp"

// Macro to define a test case
#define TEST_CASE(test_name) void test_name()

// Macro to check a condition and report success or failure
#define CHECK(condition) \
    do { \
        if (!(condition)) { \
            std::cerr << "TEST FAILED: " << #condition << " in " << __FUNCTION__ << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            return; \
        } \
    } while (0)

// Function to run a test case
void run_test(void (*test)(), const std::string& test_name) {
    std::cout << "Running TEST " << test_name << "..." << std::endl;
    test();
    std::cout  <<"Test"<< " ended: " << test_name  << std::endl;
}

void test_taujoin_single(string kmer){
    vector<string> retvec;
    tau_join( kmer, retvec);

    cout<<"[delta join length] for k-mer: "<<kmer<<" k="<<kmer.size()<<endl;
    for(auto s : retvec){
        cout<< 2*kmer.size() - s.length()<<" "<< s<< " " <<s.length() <<endl;
    }
}
TEST_CASE(test_taujoin){
    cout<<"### If length of overlap is delta, then the extra length added is k-delta. ###"<<endl;
    // test_taujoin_single("ACACA");
    // test_taujoin_single("AAAAAAAA");
    // test_taujoin_single("ACCCCCCC");
    // test_taujoin_single("AGGA");
    // test_taujoin_single("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAG");

    // test_taujoin_single("AGCTTAAAGTAATTATCTAGGTGTCTGTATTTG");
    test_taujoin_single("AAA");

    


}

TEST_CASE(test_addition) {
    int a = 2;
    int b = 3;
    int result = a + b;
    CHECK(result == 5);
}

TEST_CASE(test_variance) {
    Thm11 thm11;
    thm11.init();
    cout<<"Fixed tau "<< thm11.fixed_tau <<endl;
    double estimate = thm11.compute_sd_estimate();

    // test against the formula
    vector<int> num_replicates_vector;
    for (int i = 1; i <= 5; ++i) {
        num_replicates_vector.push_back(i * 200);
    }

    for(int num_replicates : num_replicates_vector){
        //int num_replicates = 20;
        ResultAggregator res;
        for(int i=0; i<num_replicates; i++){
            std::string mutatedString = generateMutatedString(thm11.fixed_string_s, thm11.r);
            res.test_vals_xtau.push_back(thm11.sum_xi_tau(mutatedString, thm11.fixed_tau));
        }
        res.calculateMeanSD(res.test_vals_xtau, res.test_vals_xtau_mean, res.test_vals_xtau_sd);
        res.diffReport(num_replicates, estimate, res.test_vals_xtau_mean, res.test_vals_xtau_sd, thm11.r, thm11.L, thm11.k, thm11.fixed_tau);
    }
    
    CHECK(true);
}

TEST_CASE(test_variance_antonio) {
    VarianceTest vt(0);
    vt.doIt();
}


class Tests{
    void static TauJoin();

    /// @brief * try with a random sequence but with mutation rate of 0.
    void static RandomSeqR0();

    /**
     * @brief try to verify that the variance formula is correct. 
     * fix a tau in your string, and compare the sample variance of xtau to the
     * formula in our lemma. As the number of trials grow, the two should be 
     * approaching equality.
     */ 
    void static VarianceFormula();

    // * try to track each of the values you are adding up separately. What
    // does each of the summation give you?
    // * try to output the computed values separately for each tau.
};

int main (int argc, char* argv[]){

    // run_test(test_variance, "test_variance");
    run_test(test_taujoin, "test_taujoin");
    //run_test(test_variance_antonio, "VarianceTestAntonio");



    return 0;
}

