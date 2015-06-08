
/*
ProposalEngine.hpp - This file is part of the Coalescenator (v1.0.0)


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


#ifndef __coalescenator__ProposalEngine_hpp
#define __coalescenator__ProposalEngine_hpp

#include <random>

#include "Eigen/Dense"

#include "TypeContainer.hpp"
#include "Utils.hpp"
#include "Matrix3d.hpp"
#include "SamplerOutput.hpp"
#include "LikelihoodSamples.hpp"
#include "ConditionalSamplingHMM.hpp"

using namespace std;

class ProposalEngine {

	public:

	   ProposalEngine(mt19937 *, ConditionalSamplingHMM *);
	   uint calculateLikelihood(SamplerOutput<LikelihoodSamplesType> *, ModelParameters, vector<AncestralTypeMaps>, vector<double>, vector<double>, pair<uint,uint>) ;

	private:

        mt19937 * mt_prng;
        ConditionalSamplingHMM * csHMM;

        uniform_real_distribution<double> sample_uniform_01;

    	double sampleEventTime(double);
    	double fEventTime(double, double);
        double probabilityNoEvent(double, double);
        double probabilityCombinatorics(TypeContainer &, TypeContainer &);
        double probabilitySamplePermutation(AncestralTypeMaps &);
        double piSMCJointly(Haplotype, Haplotype, TypeContainer &, uint);
        double sample_backwardMove_forwardTransition(double , double , double , double , TypeContainer &,
                                            Matrix3d<double> *, vector<double> &, vector<double> &, vector<double> &,
                                            Haplotype , double , uint, double, bool &, bool afterLastCollection);

        uint sampleFromVector(vector<double> Vector);
        double vectorSum(vector<double> Vector);

        Haplotype combineHaplotypes_AandB(Haplotype, Haplotype);
        Haplotype changeHaplotypeCtoA (Haplotype, uint);
        Haplotype changeHaplotypeCtoB (Haplotype, uint);
        pair <Haplotype,Haplotype> splitHaplotypeC (Haplotype, uint);

        double conditionalPiHat(uint, double, double, double, TypeContainer &, vector<uint> &);

        double piSMC(Haplotype, TypeContainer &, uint);

};

#endif
