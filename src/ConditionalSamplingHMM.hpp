
/*
ConditionalSamplingHMM.hpp - This file is part of the Coalescenator (v1.0.0)


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


#ifndef __coalescenator__ConditionalSamplingHMM_hpp
#define __coalescenator__ConditionalSamplingHMM_hpp

#include <vector>
#include <unordered_map>

#include <Eigen/Dense>

#include "Utils.hpp"
#include "TypeContainer.hpp"

using namespace std;

class ConditionalSamplingHMM {

	private:

		vector <double> quadratures;
		vector <double> weights;
		vector <double> r_rates;
		vector <double> thetas;

		struct HMMConstants_transition {

			vector<double> y;
			Eigen::ArrayXXd z;
		};

		Eigen::ArrayXXd e_binomial;

		unordered_map < vector<uint>, HMMConstants_transition, VectorHash<uint> > constant_container_transition;
		unordered_map < vector<uint>, Eigen::ArrayXXd, VectorHash<uint> > constant_container_emission;

		void calcConstantsTransition(uint, uint);
		void calcConstantsEmission(uint, uint, uint);

		double initialDensity(uint, uint, double);
		double calculatePartialHaplotypeEmissionProbability(uint, vector<uint> &, vector<int> &, vector<int> &);

	public:

		ConditionalSamplingHMM(uint, vector<double> &, vector<double> &, uint);
		double piSMC(Haplotype &, TypeContainer &, uint);
};

#endif
