
/*
SamplerOutput.cpp - This file is part of the Coalescenator (v1.0.0)


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

#include "SamplerOutput.hpp"
#include "Utils.hpp"
#include "Matrix3d.hpp" 
#include "LikelihoodSamples.hpp"

using namespace std;


template <typename ValueType>
SamplerOutput<ValueType>::SamplerOutput(const uint num_mus_in, const uint num_rhos_in, const uint num_lambdas_in, const uint num_samples_in, const uint num_collection_times_in) : num_mus(num_mus_in), num_rhos(num_rhos_in), num_lambdas(num_lambdas_in), num_samples(num_samples_in), num_collection_times(num_collection_times_in), likelihood_(Matrix3d<ValueType>(num_mus, num_rhos, num_lambdas, num_samples)), likelihoodSq_(Matrix3d<ValueType>(num_mus, num_rhos, num_lambdas, num_samples)), tmrca_(Matrix3d<ValueType>(num_mus, num_rhos, num_lambdas, num_samples)), lineages_remaining_(vector<Matrix3d<ValueType> >(num_collection_times, Matrix3d<ValueType>(num_mus, num_rhos, num_lambdas, num_samples))) {

}


template <typename ValueType>
Matrix3d<ValueType> & SamplerOutput<ValueType>::likelihood() {

	return likelihood_;
}


template <typename ValueType>
Matrix3d<ValueType> & SamplerOutput<ValueType>::likelihoodSq() {
    
	return likelihoodSq_;
}


template <typename ValueType>
Matrix3d<ValueType> & SamplerOutput<ValueType>::tmrca() {

	return tmrca_;
}


template <typename ValueType>
vector<Matrix3d<ValueType> > & SamplerOutput<ValueType>::lineages_remaining() {

	return lineages_remaining_;
}


template <>
void SamplerOutput<LikelihoodSamplesType>::addSample(Matrix3d<double> & likelihood_matrix_in, double waitingTimeSum, vector<double> & lineages_left) {

    likelihood_.addMatrixElementWise(likelihood_matrix_in); // linear add ok as log-likelihoods are only cumulated here, they are summed in the log domain later
    likelihoodSq_.multiplyValueAndAddMatrixElementWise(likelihood_matrix_in, 2);

    //this needed corection anyway for the rightmost mrca term

    assert(waitingTimeSum > 0);

    tmrca_.addAndAddMatrixElementWise(likelihood_matrix_in, log(waitingTimeSum));
    assert(lineages_remaining_.size() == lineages_left.size());

    for (uint l = 0; l < lineages_left.size(); l++) {

        assert(lineages_left[l] > 0);
        lineages_remaining_[l].addAndAddMatrixElementWise(likelihood_matrix_in, log(lineages_left[l]));
    }
}


template <>
void SamplerOutput<LikelihoodSamplesType>::addSamples(SamplerOutput<LikelihoodSamplesType> & samples) {

    likelihood_.addMatrixElementWise(samples.likelihood());
    likelihoodSq_.addMatrixElementWise(samples.likelihoodSq());

    //this needed corection anyway for the rightmost mrca term
    tmrca_.addMatrixElementWise(samples.tmrca());

    assert(lineages_remaining_.size() == samples.lineages_remaining_.size());

    for (uint l = 0; l < samples.lineages_remaining_.size(); l++) {

        lineages_remaining_[l].addMatrixElementWise(samples.lineages_remaining()[l]);
    }
}


template <>
void SamplerOutput<double>::sumSamples(SamplerOutput<LikelihoodSamplesType> & all_samples) {

    likelihood_.addSummedMatrixElementWise(all_samples.likelihood());
    likelihoodSq_.addSummedMatrixElementWise(all_samples.likelihoodSq());
    tmrca_.addSummedMatrixElementWise(all_samples.tmrca());
    tmrca_.subtractMatrixElementWise(likelihood_);

    assert(num_collection_times == all_samples.lineages_remaining().size());
    assert(num_collection_times == lineages_remaining_.size());

    for (uint l = 0; l < num_collection_times; l++) {

        lineages_remaining_[l].addSummedMatrixElementWise(all_samples.lineages_remaining()[l]);
        lineages_remaining_[l].subtractMatrixElementWise(likelihood_);
    }

    likelihood_.subtractConstant(num_samples);
    likelihoodSq_.subtractConstant(num_samples);
}


template class SamplerOutput<double>;
template class SamplerOutput<LikelihoodSamplesType>;
