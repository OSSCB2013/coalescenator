
/*
LikelihoodSamples.cpp - This file is part of the Coalescenator (v1.0.0)


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


#include <vector>
#include <assert.h>
#include <set>

#include "LikelihoodSamples.hpp"
#include "Utils.hpp"

using namespace std;

LikelihoodSamples::LikelihoodSamples(const uint num_samples_in) : num_samples(num_samples_in) {}

void LikelihoodSamples::operator += (const double & value) {
    
    assert(std::isfinite(value));
    likelihood_samples_.insert(value);
}


void LikelihoodSamples::operator += (const LikelihoodSamples & other) {

    likelihood_samples_.insert(other.likelihood_samples_.begin(), other.likelihood_samples_.end());
}


LikelihoodSamplesSum::LikelihoodSamplesSum(const uint num_samples_in) : LikelihoodSamples(num_samples_in) {}

double LikelihoodSamplesSum::sumSamples() {

    assert(num_samples == likelihood_samples_.size());

    auto sit = likelihood_samples_.rbegin();
    double sample_sum = *sit;
    sit++;

    while (sit != likelihood_samples_.rend()) {
        sample_sum = logAddition(sample_sum, *sit);
        sit++;
    }

    return sample_sum;
}


LikelihoodSamplesKahan::LikelihoodSamplesKahan(const uint num_samples_in) : LikelihoodSamples(num_samples_in) {}

double LikelihoodSamplesKahan::sumSamples() {

    auto sit = likelihood_samples_.begin();
    double sample_sum = *sit;
    sit++;

    double compensation = -pow(10,(double)300);
    assert(compensation + pow(10,(double)20) < sample_sum);

    while (sit != likelihood_samples_.end()) {

        double compensated_value = logSubtraction(*sit, compensation);
        double temporary_sum = logAddition(sample_sum, compensated_value);
        double sum_difference = logSubtraction(temporary_sum, sample_sum);
        compensation = logSubtraction(sum_difference, compensated_value);
        
        sample_sum = temporary_sum;
        sit++;
    }

    return sample_sum;
}