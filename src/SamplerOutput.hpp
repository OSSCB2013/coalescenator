
/*
SamplerOutput.hpp - This file is part of the Coalescenator (v1.0.0)


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


#ifndef __coalescenator__SamplerOutput_hpp
#define __coalescenator__SamplerOutput_hpp

#include <vector>

#include "Utils.hpp"
#include "Matrix3d.hpp" 
#include "LikelihoodSamples.hpp"

using namespace std;

template <typename ValueType>
class SamplerOutput {

    public: 

    	SamplerOutput(const uint, const uint, const uint, const uint, const uint);

    	Matrix3d<ValueType> & likelihood();
        Matrix3d<ValueType> & likelihoodSq();
    	Matrix3d<ValueType> & tmrca();
    	vector<Matrix3d<ValueType> > & lineages_remaining();

        void addSample(Matrix3d<double> &, double, vector<double> &);
        void addSamples(SamplerOutput<LikelihoodSamplesType> &);
        void sumSamples(SamplerOutput<LikelihoodSamplesType> &);

	private: 
    	
    	const uint num_mus;
    	const uint num_rhos;
    	const uint num_lambdas;
        const uint num_samples;
        const uint num_collection_times;

    	Matrix3d<ValueType> likelihood_;
        Matrix3d<ValueType> likelihoodSq_;
    	Matrix3d<ValueType> tmrca_;
    	vector<Matrix3d<ValueType> > lineages_remaining_;
};


#endif