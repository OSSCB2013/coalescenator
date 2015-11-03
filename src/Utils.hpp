
/*
Utils.hpp - This file is part of the Coalescenator (v1.0.0)


The MIT License (MIT)

Copyright (c) 2015 Kevin Dialdestoro, Jonas Andreas Sibbesen, Lasse Maretty and Paul Jenkins

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/


#ifndef __coalescenator__Utils_hpp
#define __coalescenator__Utils_hpp

#include <string>
#include <vector>
#include <time.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <unordered_map>
#include <limits>
#include <functional>

#include "boost/functional/hash.hpp"

using namespace std;

const double double_underflow = numeric_limits<double>::min();
const double double_overflow = numeric_limits<double>::max();
const double double_precision = numeric_limits<double>::epsilon();

typedef unsigned int uint;
typedef unsigned char uchar;

template <typename T>
size_t hash_value(T val) {

    return std::hash<T>(val);
}

template <typename T>
struct VectorHash {

    size_t operator()(const vector<T> & vec) const {

      size_t seed = 0;
      hash<T> hash_fn;

      // From boost 1.58.0 - functional/hash/hash.hpps: void hash_combine_impl(SizeT& seed, SizeT value)
      for(auto elem : vec) {

          seed ^= hash_fn(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      }

      return seed;
    }
};

// template <typename T>
// struct VectorHash {
//
//     size_t operator()(const vector<T> & c) const {
//
//         return boost::hash_range(c.begin(), c.end());
//     }
// };


inline double logAddition(const double log_summand1, const double log_summand2) {

    assert(std::isfinite(log_summand1));
    assert(std::isfinite(log_summand2));

    if (log_summand1 < log_summand2) {

        double log_sum = log_summand2 + log(1 + exp(log_summand1 - log_summand2));
        assert(std::isfinite(log_sum));

        return log_sum;

    } else {

        double log_sum = log_summand1 + log(1 + exp(log_summand2 - log_summand1));
        assert(std::isfinite(log_sum));

        return log_sum;
    }
}


inline double logSubtraction(const double log_summand1, const double log_summand2) {

    assert(std::isfinite(log_summand1));
    assert(std::isfinite(log_summand2));

    double log_sum = log_summand1 + log(1 - exp(log_summand2 - log_summand1));

    if (!std::isfinite(log_sum)) {

        log_sum = -pow(10,(double)300);
    }

    return log_sum;
}


typedef unordered_map<vector<bool>, uint, VectorHash<bool> > TypeMap;

// Container for option variables (all parameters used in the code)
struct OptionsContainer {

    string fasta_list;
    string time_list;
    string viral_list;
    double generation_time;

    string coordinate_list;
    string output_file_prefix;
    uint loci_size;

    uint num_threads;
    uint seed;

    string scale_list;
    string mutation_list;
    string rho_list;
    double scale_driving;
    double mutation_driving;
    double rho_driving;

    uint mc_samples;
    uint qadrature_points;

    uint min_dist;
    uint max_dist;
};

struct ModelParameters {

    double driving_mu;
    double driving_lambda;
    double driving_rho;
    vector<double> target_mu;
    vector<double> target_lambda;
    vector<double> target_rho;
    double tau;
};

struct Haplotype {

    vector<bool> type;
    string gamma;
};

struct AncestralTypeMaps {

    TypeMap A;
    TypeMap B;
    TypeMap C;
};

// Vector flush capability
template < typename T >
std::ostream& operator << (std::ostream& os, typename std::vector<T>& v) {

    for (typename vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii) {

		os << scientific << setprecision(6)  << *ii << "  ";
	}

    return os;
}

// 2D-vector flush capability
template < typename T >
std::ostream& operator << (std::ostream& os, typename std::vector<std::vector<T> >& v) {

    for (typename vector<vector<T> >:: iterator ii = v.begin(); ii != v.end(); ++ii) {

        os << *ii << "\n";
    }

    return os;
}


inline bool double_compare(double a, double b) {

    return ((a == b) or (abs(a - b) < abs(min(a, b)) * double_precision));
}

inline string getLocalTime () {

    time_t now = time(NULL);
    struct tm * lt;
    lt = localtime (&now);

    stringstream output_string;

    if (lt->tm_mday < 10) {

        output_string << "0" << lt->tm_mday << "/";

    } else {

        output_string << lt->tm_mday << "/";
    }

    if ((1 + lt->tm_mon) < 10) {

        output_string << "0" << (1 + lt->tm_mon) << "/" << (1900 + lt->tm_year) << " ";

    } else {

        output_string << (1 + lt->tm_mon) << "/" << (1900 + lt->tm_year) << " ";
    }

    if (lt->tm_hour < 10) {

        output_string << "0" << lt->tm_hour << ":";

    } else {

        output_string << lt->tm_hour << ":";
    }

    if (lt->tm_min < 10) {

        output_string << "0" << lt->tm_min << ":";

    } else {

        output_string << lt->tm_min << ":";
    }

    if (lt->tm_sec < 10) {

        output_string << "0" << lt->tm_sec;

    } else {

        output_string << lt->tm_sec;
    }

    return output_string.str();
}

#endif
