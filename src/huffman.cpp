#include <vector>
#include <vector>
#include <iostream>
#include <cassert>

using u32 = std::uint32_t;

std::vector<u32> construct_canonical_code( std::vector<u32> const& lengths ){
    unsigned int size = lengths.size();
    std::vector< unsigned int > length_counts(16,0); //Lengths must be less than 16 for DEFLATE
    u32 max_length = 0;
    for(auto i: lengths){
        assert(i <= 15);
        length_counts.at(i)++;
        max_length = std::max(i, max_length);
    }
    length_counts[0] = 0; //Disregard any codes with alleged zero length

    std::vector< u32 > result_codes(size,0);

    //The algorithm below follows the pseudocode in RFC 1951
    std::vector< unsigned int > next_code(size,0);
    {
        //Step 1: Determine the first code for each length
        unsigned int code = 0;
        for(unsigned int i = 1; i <= max_length; i++){
            code = (code+length_counts.at(i-1))<<1;
            next_code.at(i) = code;
        }
    }
    {
        //Step 2: Assign the code for each symbol, with codes of the same length being
        //        consecutive and ordered lexicographically by the symbol to which they are assigned.
        for(unsigned int symbol = 0; symbol < size; symbol++){
            unsigned int length = lengths.at(symbol);
            if (length > 0)
                result_codes.at(symbol) = next_code.at(length)++;
        }  
    } 
    return result_codes;
}

static bool sortByProb(const std::pair<std::vector<int>, double>& a, const std::pair <std::vector<int>, double>& b){
    return (a.second < b.second);
}

void printPairVector(std::vector < std::pair <std::vector<int>, double> >& v){
    std::cerr << v.size()<< std::endl;
    for (const auto& pair : v) {
        std::cerr << "symbol: ";
        for (const auto& symbol : pair.first) {
            std::cerr << symbol << " ";
        }
        std::cerr << "prob: " << pair.second << std::endl;
    }
}

void printLengths(const std::vector<int>& symbols, const std::vector<int>& lengths_table){
    std::cerr << "lengths table: " << std::endl;
    for(u32 idx = 0; idx < symbols.size(); idx++){
        std::cerr << symbols.at(idx) << "    " << lengths_table.at(idx) << std::endl;
    }
}

template<typename T>
void appendVectors(std::vector<T>& dest, const std::vector<T>& src1, const std::vector<T>& src2) {
    dest.insert(dest.end(), src1.begin(), src1.end());
    dest.insert(dest.end(), src2.begin(), src2.end());
}

u32 countOccurances(std::vector < std::pair <std::vector<int>, double> >& probabilities, u32 symbol){
    u32 count {0};
    for (const auto& pair : probabilities){         // traverse properties vector
        for (const auto& curr_symbol : pair.first) {     // traverse set of symbol
            if (curr_symbol == symbol){
                ++count;
            }
        }
    }
    return count;
}

std::vector<int> package_merge(const std::vector <std::pair <std::vector<int>, double>> probabilities, u32 num_symbols){
    u32 size = probabilities.size();
    // get symbols
    std::vector<int> symbols;
    for(auto pair : probabilities){
        symbols.push_back(pair.first.at(0));
    }
    std::vector <std::pair <std::vector<int>, double>> current {probabilities};

    while(current.size() < (2*size)-2){
        // sort with increasing probabilities
        sort(current.begin(), current.end(), sortByProb);
        std::cerr << "---- packages ----" << std::endl;
        printPairVector(current);
        // if odd number of packages -> discard last (if it not one of the original packages)
        if(current.size()%2 != 0){
            current.pop_back();
        }

        std::vector < std::pair <std::vector<int>, double> > packages {};
        for(u32 idx = 0; idx < current.size()-1; idx += 2){
            // create package
            std::vector<int> v {};
            appendVectors(v, current.at(idx).first, current.at(idx+1).first);

            double prob = current.at(idx).second + current.at(idx+1).second;

            // merge
            packages.push_back({v, prob});
        }

        current.clear();
        appendVectors(current, probabilities, packages);
    }

    std::vector<int> lengths_table {};
    for(int symbol : symbols){
        lengths_table.push_back(countOccurances(current, symbol));
    }

    printLengths(symbols, lengths_table);

    return lengths_table;
}

bool are_valid_lengths(std::vector<u32> lengths){
    double sum = 0;
    for(u32 l : lengths){
        sum = sum + (1/(1 << l));
    }
    return sum <= 1;
}


int main()
{
 
    std::vector <std::pair <std::vector<int>, double>> probabilities {
        {{-100},    0.004},
        {{-5},      0.002},
        {{-4},      0.004},
        {{-3},      0.010},
        {{-2},      0.033},
        {{-1},      0.194},
        {{0},       0.245},
        {{1},       0.191},
        {{2},       0.034},
        {{3},       0.010},
        {{4},       0.004},
        {{5},       0.002},
        {{100},     0.004},
        {{120},     0.018},
        {{150},     0.245}
    };

    std::vector <std::pair <std::vector<int>, double>> probabilities2 {
        {{-1},  0.11},
        {{-2},  0.15},
        {{-3},  0.45},
        {{-4},  0.08},
        {{-5},  0.01},
        {{-6},  0.20}
    };

    std::vector<u32> lenghts = {4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3};
    std::vector<u32> option1 = {9, 9, 8, 7, 5, 2, 2, 3, 6, 7, 8, 9, 9, 5, 2};

    std::cerr << are_valid_lengths(option1) <<std::endl;
    std::vector<u32> encodings = construct_canonical_code(option1);

    for(u32 idx = 0; idx < encodings.size(); idx++){
        std::cerr << encodings.at(idx) << std::endl;
    }

    // package_merge(probabilities2, probabilities2.size());
    // package_merge(probabilities, probabilities.size());
 
    return 0;
}
