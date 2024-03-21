#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>
#include<sstream>
#include <fstream>
#include <vector>
#include <bitset>

using namespace std;

//vector<string> kmers;
vector<int> counts;
vector<uint64_t> kmer_binary_values;

int K = 30;
string zodd="";
string zeven="";
uint64_t zodd_t, zeven_t;

/// @brief murmurhash3 64-bit finalizer
/// @param key 
/// @return 
inline uint64_t hash64(uint64_t key )
{
    int k = K;
    uint64_t mask = (1ULL<<k*2) - 1;
    key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
    key = key ^ key >> 14;
    key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
    key = key ^ key >> 28;
    key = (key + (key << 31)) & mask;
    return key;
}


/*
/// @brief returns 0 to 1 real value for a k-mer
/// @param kmer 
/// @return 
double kmer_to_real_value(const std::string& kmer) {
    int k = kmer.length();
    double binary_value = 0;
    for (int i = 0; i < k; ++i) {
        switch (kmer[i]) {
            case 'A':
                binary_value = binary_value * 4; // 'A' is represented as 00 in binary
                break;
            case 'C':
                binary_value = binary_value * 4 + 1; // 'C' is represented as 01 in binary
                break;
            case 'G':
                binary_value = binary_value * 4 + 2; // 'G' is represented as 10 in binary
                break;
            case 'T':
                binary_value = binary_value * 4 + 3; // 'T' is represented as 11 in binary
                break;
            default:
                throw std::invalid_argument("Invalid character in k-mer: " + kmer[i]);
        }
    }

    // Normalize binary_value to range [0, 1]
    //K = k;
    double real_value = hash64(binary_value) / std::pow(4, k);

    return real_value;
}*/

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

uint64_t stringToBinary(string kmer){
    string bv_line = encode_kmer(kmer);
    //int K = 30;
    uint64_t curr_bv_lo = std::stoull(bv_line.substr(0,std::min(64, 2*K)), nullptr, 2);
    // uint64_t curr_bv_hi = 0;
    // if(K >= 64){
    //     curr_bv_hi = std::stoull(bv_line.substr(64,bv_line.length()-64), nullptr, 2);
    // } 
    //int hd = hammingDistance(prev_bv_hi, curr_bv_hi);
    cout<<bv_line<<" "<<curr_bv_lo<<endl;

    return curr_bv_lo;
} 

double kmer_to_real_value(const std::string& kmer) {
    return hash64(stringToBinary(kmer))*1.0/std::pow(4, K);
}

double kmer_binary_to_real_value(uint64_t kmer_binary) {
    return hash64(kmer_binary)*1.0/std::pow(4, K);
}

//write a function that reads all lines from a file and put in a vector of strings  
void read_kmers(string file_name){
    ifstream ifs(file_name);
    string line;
    // vector<string> lines;
    while (getline(ifs, line))
    {
        stringstream ss(line);
        string kmer;
        int count;

        if (ss >> kmer >> count) {
            if(kmer_binary_values.empty()){
                K = kmer.length();
            }

            //kmers.push_back(kmer);
            kmer_binary_values.push_back(stringToBinary(kmer));
            counts.push_back(count);
        }

    }
    
    ifs.close();
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

void getSketch(double threshold, vector<uint64_t>& sketch_kmers, vector<int>& sketch_counts){//threshold is the minimum value of the sketch to be considered
    for(int i=0; i< kmer_binary_values.size(); i++){
        uint64_t kmer = kmer_binary_values[i];
        if(kmer_binary_to_real_value(kmer) <= threshold){
            //cout << kmer << " "<< kmer_to_real_value(kmer)<< endl;
            sketch_kmers.push_back(kmer_binary_values[i]);
            sketch_counts.push_back(counts[i]);
        }  
    }
    cout<< "Number of original kmers: " << kmer_binary_values.size() << endl;
    cout<< "Number of sketch kmers: " << sketch_kmers.size() << endl;
}


int hammingDistance (uint64_t x, uint64_t y) {

// Let x and y be the 2-bit encoding of two kmers.

// let zeven = (x with every even position zeroed out) XOR (y with every
// even position zeroed out)
// let zodd be similarly defined
// HD(x,y) = popcount(zeven OR zodd)

    // std::bitset<64> xo(x&zodd_t); // Assuming a 64-bit binary representation
    // std::bitset<64> yo(y&zodd_t); // Assuming a 64-bit binary representation
    // std::bitset<64> xb(x); // Assuming a 64-bit binary representation
    // std::bitset<64> yb(y); // Assuming a 64-bit binary representation

    //     //cout<<xb.to_string()<<" "<<yb.to_string()<<endl;
        uint64_t oddt = (x&zodd_t)^(y&zodd_t);
        uint64_t event = (x&zeven_t)^(y&zeven_t);

    // std::bitset<64> oddtt(oddt); // Assuming a 64-bit binary representation

    //     cout<<xb.to_string()<<" "<<"x"<<endl;
    //     cout<<yb.to_string()<<" "<<"y"<<endl;

    //     cout<<xo.to_string()<<" "<<"x&zo"<<endl;
    //     cout<<yo.to_string()<<" "<<"y&zo"<<endl;

    //     cout<<oddtt.to_string()<<" "<<"oddxor"<<endl;
    //     cout<<oddtt.to_string()<<" "<<"oddxor"<<endl;

        // 0101010101010101010101010101010101010101010101010101010101010101

		// uint64_t res = x ^ y;
		// return __builtin_popcountll (res) ;
        
		return __builtin_popcountll (oddt | (event>>1)) ;
	}




inline unsigned int num_pair(int n){
    return (n*n - n)/2;
}

void getHDHistogram(double threshold){
    vector<uint64_t> sketch_kmers_binary;
    vector<int> sketch_kmers_count;

    getSketch(threshold, sketch_kmers_binary, sketch_kmers_count); 
    vector<int> hd_histogram(K+1, 0);
    for (int i = 0; i < sketch_kmers_binary.size(); i++){
        for (int j = 0; j < sketch_kmers_binary.size(); j++){
            int hd = hammingDistance(sketch_kmers_binary[i], sketch_kmers_binary[j]);
            hd_histogram[hd] += sketch_kmers_count[j];
        }
        // if(counts[i] > 1){
        //     hd_histogram[0] += num_pair(counts[i]); 
        // }
        // for (int j = i+1; j < sketch_kmers_binary.size(); j++){
        //     //int hd = hammingDistance(stringToBinary(kmers[i]), stringToBinary(kmers[j]));
        //     //int hd = hammingDist(kmers[i], kmers[j]);
        //     int hd = hammingDistance(sketch_kmers_binary[i], sketch_kmers_binary[j]);
        //     //hd_histogram[hd] += 1; // if count 1
        //     hd_histogram[hd] += counts[i]*counts[j];
        // }
    }
    for (int i = 0; i < hd_histogram.size(); i++){
        cout << i << " " << hd_histogram[i] << endl;
    }

}


int main(int argc, char *argv[]) {

    // cout<<kmer_to_real_value("ATTTTTAAAAAAAATATATATGGATATATA");
    // exit(1);
    // string kmer1="CTTTTTAAAAAAAATATATATGGATATATA";
    // string kmer2="CTTTTTGAAAAAAATATATATGGATATAAT";
    // cout<<hammingDistance(stringToBinary(kmer1), stringToBinary(kmer2));
    // exit(1);
    read_kmers("kmers.txt"); // K IS SET


    for (int i = 0; i < K+2; ++i) { //i<32 for K=30
        zodd += "01";
    }
    for (int i = 0; i < K+2; ++i) {
        zeven += "10";
    }
    zodd_t = std::stoull(zodd, nullptr, 2);
    zeven_t = std::stoull(zeven, nullptr, 2);


    double value = std::stod(argv[1]);
    getHDHistogram(value);

    //std::string kmer = "AAAAAAAAATAAAAAAAAAA";
    //double real_value = kmer_to_real_value(kmer);
    //std::cout << "Real value for k-mer " << kmer << ": " << real_value << std::endl;

    return 0;
}
