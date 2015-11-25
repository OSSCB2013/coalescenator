
/*
ProposalEngine.cpp - This file is part of the Coalescenator (v1.0.0)


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


#include <iostream>
#include <vector>
#include <math.h>

#include <boost/random/exponential_distribution.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "ProposalEngine.hpp"
#include "TypeContainer.hpp"
#include "Utils.hpp"
#include "ConditionalSamplingHMM.hpp"
#include "LikelihoodSamples.hpp"

using namespace std;

ProposalEngine::ProposalEngine(mt19937 * mt_prng_in, ConditionalSamplingHMM * csHMM_in) {

    mt_prng = mt_prng_in;
    csHMM = csHMM_in;
    sample_uniform_01 = uniform_real_distribution<double>(0,1);
}

double ProposalEngine::sampleEventTime(double rate) {

    double uniform_sample = sample_uniform_01(*mt_prng);

    return -log(uniform_sample) / rate;

}

double ProposalEngine::fEventTime(double t, double rate) {

    return log(rate) - rate*t;

}

double ProposalEngine::probabilityNoEvent(double t, double rate) {

    return -rate*t;

}

double ProposalEngine::probabilityCombinatorics(TypeContainer &H, TypeContainer &H_beforeCollection) {

    vector < pair <uint, uint> > type_map_overlap = H.getTypeMapOverlap();
    assert (type_map_overlap.size() > 0);

    assert (type_map_overlap.size() <= H.getUniqueCount());
    
    uint lineage_count_H_A = H.getLineageCountA();
    double combi_factor = -log(boost::math::binomial_coefficient<double>(lineage_count_H_A, lineage_count_H_A - H_beforeCollection.getLineageCountA()));

    uint lineage_count_H_B = H.getLineageCountB();
    combi_factor -= log(boost::math::binomial_coefficient<double>(lineage_count_H_B, lineage_count_H_B - H_beforeCollection.getLineageCountB()));

    uint lineage_count_H_C = H.getLineageCountC();
    combi_factor -= log(boost::math::binomial_coefficient<double>(lineage_count_H_C, lineage_count_H_C - H_beforeCollection.getLineageCountC()));

    for (uint i = 0; i < type_map_overlap.size(); i++) {

        combi_factor += log(boost::math::binomial_coefficient<double>(type_map_overlap[i].first, type_map_overlap[i].second));
    }

    return combi_factor;

}

double ProposalEngine::probabilitySamplePermutation(AncestralTypeMaps & type_maps) {

    double probability = 0;

    uint total_A_count = 0;
    for (TypeMap::iterator it = type_maps.A.begin(); it != type_maps.A.end(); it++){

        assert (it->second > 0);
        total_A_count += it->second;
        probability += boost::math::lgamma(it->second + 1);
    }

    probability -= boost::math::lgamma(total_A_count + 1);

    uint total_B_count = 0;
    for (TypeMap::iterator it = type_maps.B.begin(); it != type_maps.B.end(); it++){

        assert (it->second > 0);
        total_B_count += it->second;
        probability += boost::math::lgamma(it->second + 1);
    }

    probability -= boost::math::lgamma(total_B_count + 1);

    uint total_C_count = 0;
    for (TypeMap::iterator it = type_maps.C.begin(); it != type_maps.C.end(); it++){

        assert (it->second > 0);
        total_C_count += it->second;
        probability += boost::math::lgamma(it->second + 1);
    }

    probability -= boost::math::lgamma(total_C_count + 1);

    assert(probability <= 0);

    return probability;
}

double ProposalEngine::piSMCJointly(Haplotype HaplotypeA, Haplotype HaplotypeB, TypeContainer &H, uint c) {

    double piSMC_temp = csHMM->piSMC(HaplotypeA, H, c);
    H.addType(HaplotypeA);
    double first_piSMC = 0.5 * piSMC_temp * csHMM->piSMC(HaplotypeB, H, c);
    H.removeType(HaplotypeA);

    piSMC_temp = csHMM->piSMC(HaplotypeB, H, c);
    H.addType(HaplotypeB);
    double second_piSMC = 0.5 * piSMC_temp * csHMM->piSMC(HaplotypeA, H, c);
    H.removeType(HaplotypeB);

    return first_piSMC + second_piSMC;

}



double ProposalEngine::sample_backwardMove_forwardTransition(double driving_N, double driving_theta, double driving_rho, double tau, TypeContainer &H,
                                            Matrix3d<double> * likelihood_matrix,  vector<double> &target_mu, vector<double> &target_r, vector<double> &target_lambda,
                                            Haplotype chosenHaplotype, double waitingTime, uint c, double V, bool &postCollection, bool afterLastCollection)
{
    uint totalSequenceLength = H.getTotalSequenceLength();
    uint locus_len_A = H.getLocusLengthA();
    uint locus_len_B = H.getLocusLengthB();

    double driving_thetaA = locus_len_A * driving_theta;
    double driving_thetaB = locus_len_B * driving_theta;
    double driving_thetaAB = totalSequenceLength * driving_theta;

    Haplotype chosenHaplotype2;
    Haplotype updatedHaplotype;
    pair<Haplotype, Haplotype> updatedHaplotype_Recombination;

    uint Cij = H.getTypeOccurrenceC(chosenHaplotype);
    uint Ci = H.getTypeOccurrenceCa(chosenHaplotype);
    uint Cj = H.getTypeOccurrenceCb(chosenHaplotype);
    uint Ai = H.getTypeOccurrenceA(chosenHaplotype);
    uint Bj = H.getTypeOccurrenceB(chosenHaplotype);

    uint n = H.getLineageCount();
    uint C = H.getLineageCountC();
    uint A = H.getLineageCountA();
    uint B = H.getLineageCountB();

    assert (n == A + B + C);

    double event_constant; //denote one of the 8 events constant //CoalescenceC,CoalescenceA,CoalescenceB,CoalescenceAB,Mutation,MutationA,MutationB,Recombination
    string event_string;

    vector<double> probabilityVector;
    uint sampledIndex;

    H.removeType(chosenHaplotype);
    double no_event_pi = csHMM->piSMC(chosenHaplotype, H, c);

    if (chosenHaplotype.gamma == "C")  //a sequence of the form (i,j) is chosen
    {
        probabilityVector.resize(totalSequenceLength+4, 0);

        assert (Cij > 0);
        assert (Ci > 0);
        assert (Cj > 0);

        //mutation
        for (uint i = 0; i < totalSequenceLength; i++) {

            assert (chosenHaplotype.type.size() == totalSequenceLength);
            updatedHaplotype = chosenHaplotype;
            updatedHaplotype.type[i] = (updatedHaplotype.type[i]+1)%2;

            double piSMCSampledMutation = csHMM->piSMC(updatedHaplotype, H, c);
            probabilityVector[i] = driving_theta*Cij*piSMCSampledMutation/no_event_pi;
            assert (probabilityVector[i] >= 0);

        }

        //coalescence of two (i,j)s
        if (Cij > 1) {

            probabilityVector[totalSequenceLength] = Cij*(Cij-1)/no_event_pi;
            assert (probabilityVector[totalSequenceLength] > 0);

        } else {

            probabilityVector[totalSequenceLength] = 0;

        }

        //coalescence of (i,j) with (i,*)
        if (Ai > 0) {

            chosenHaplotype2 = changeHaplotypeCtoA(chosenHaplotype, locus_len_A);
            assert (chosenHaplotype2.type.size() == locus_len_A);

            H.addType(chosenHaplotype);
            H.removeType(chosenHaplotype2);

            probabilityVector[totalSequenceLength+1] = Ai*Cij/csHMM->piSMC(chosenHaplotype2, H,c);
            assert (probabilityVector[totalSequenceLength + 1] > 0);
            
            H.addType(chosenHaplotype2);
            H.removeType(chosenHaplotype);


        } else {

            probabilityVector[totalSequenceLength+1] = 0;

        }

        //coalescence of (i,j) with (*,j)
        if (Bj > 0) {

            chosenHaplotype2 = changeHaplotypeCtoB(chosenHaplotype, locus_len_A);
            assert (chosenHaplotype2.type.size() == locus_len_B);
            
            H.addType(chosenHaplotype);
            H.removeType(chosenHaplotype2);

            probabilityVector[totalSequenceLength+2] = Bj*Cij/csHMM->piSMC(chosenHaplotype2, H,c);
            assert (probabilityVector[totalSequenceLength + 2] > 0);
            
            H.addType(chosenHaplotype2);
            H.removeType(chosenHaplotype);

        } else {

            probabilityVector[totalSequenceLength+2] = 0;

        }

        //recombination
        updatedHaplotype_Recombination = splitHaplotypeC(chosenHaplotype, locus_len_A);
        assert (updatedHaplotype_Recombination.first.type.size() == locus_len_A);
        assert (updatedHaplotype_Recombination.second.type.size() == locus_len_B);

        double piSMCSampling2Marginals = piSMCJointly(updatedHaplotype_Recombination.first, updatedHaplotype_Recombination.second, H, c);

        probabilityVector[totalSequenceLength+3] = driving_rho*Cij*piSMCSampling2Marginals/no_event_pi;
        assert (probabilityVector[totalSequenceLength + 3] >= 0);

        //sampling transition move
        sampledIndex = sampleFromVector(probabilityVector);

        if (sampledIndex < totalSequenceLength){ //have a mutation (i,j) to (k,j) or (i,j) to (i,l)

            updatedHaplotype = chosenHaplotype;
            updatedHaplotype.type[sampledIndex] = (updatedHaplotype.type[sampledIndex]+1) % 2;

            H.addType(updatedHaplotype);

            double Ckj = H.getTypeOccurrenceC(updatedHaplotype); //equivalent to Cij+1 when compared to Hk
            event_constant = Ckj;
            event_string = "Mutation";

        } else if (sampledIndex == totalSequenceLength){ //coalescence of two (i,j)s

            assert (Cij > 1);

            // event_constant = n*(Cij - 1);
            event_constant = C*(Cij - 1);

            event_string = "Coalescence";

        } else if (sampledIndex == totalSequenceLength + 1){ //coalescence of (i,j) with (i,*)

            assert (Ai > 0);
            assert (Ci > 0);

            chosenHaplotype2 = changeHaplotypeCtoA(chosenHaplotype, locus_len_A);
            assert (chosenHaplotype2.type.size() == locus_len_A);

            H.addType(chosenHaplotype);
            H.removeType(chosenHaplotype2);
            
            // event_constant = n*(Ai + 2*Ci - 1);
			event_constant = A*Cij;

            event_string = "Coalescence";

        } else if (sampledIndex == totalSequenceLength + 2){ //coalescence of (i,j) with (*,j)

            assert (Bj > 0);
            assert (Cj > 0);

            chosenHaplotype2 = changeHaplotypeCtoB(chosenHaplotype, locus_len_A);
            assert (chosenHaplotype2.type.size() == locus_len_B);
            
            H.addType(chosenHaplotype);
            H.removeType(chosenHaplotype2); //remove (*,j)
            
            // event_constant = n*(Bj + 2*Cj - 1);
			event_constant = B*Cij;

            event_string = "Coalescence";

        } else if (sampledIndex == totalSequenceLength + 3){ //recombination of (i,j)

            updatedHaplotype_Recombination = splitHaplotypeC(chosenHaplotype, locus_len_A);
            assert (updatedHaplotype_Recombination.first.type.size() == locus_len_A);
            assert (updatedHaplotype_Recombination.second.type.size() == locus_len_B);

            H.addType(updatedHaplotype_Recombination.first);
            H.addType(updatedHaplotype_Recombination.second);

            double Ai = H.getTypeOccurrenceA(updatedHaplotype_Recombination.first); //equivalent to Ai+1 when compared to Hk
            double Bj = H.getTypeOccurrenceB(updatedHaplotype_Recombination.second); //equivalent to Bj+1 when compared to Hk

            assert (Ai > 0);
            assert (Bj > 0);

            // event_constant = Ai*Bj/(n+1);
            event_constant = C*Ai*Bj/((A+1.0)*(B+1.0));

            event_string = "Recombination";

        } else {

            assert (1 == 0);

        }


    }

    //if sequence of the form (i,*) is chosen
    if (chosenHaplotype.gamma == "A"){

        assert (Ai > 0);
        vector<Haplotype> vectorHaplotypesB = H.getHaplotypesB();
        probabilityVector.resize(locus_len_A + vectorHaplotypesB.size() + 1, 0);

        //mutation to (i,*)
        for (uint i = 0; i < locus_len_A; i++)
        {

            assert (chosenHaplotype.type.size() == locus_len_A);
            updatedHaplotype = chosenHaplotype;
            updatedHaplotype.type[i] = (updatedHaplotype.type[i]+1)%2;
            
            double piSMCSampledMutation = csHMM->piSMC(updatedHaplotype, H,c);
            probabilityVector[i] = driving_theta*Ai*piSMCSampledMutation/no_event_pi;
            assert (probabilityVector[i] >= 0);

        }

        //coalescence (i,*) with (*,j)
        for (uint i = 0; i < vectorHaplotypesB.size(); i++){

            chosenHaplotype2 = vectorHaplotypesB[i];
            updatedHaplotype = combineHaplotypes_AandB(chosenHaplotype, chosenHaplotype2); //input (i,*), (*,j) and return (i,j)
            assert (updatedHaplotype.type.size() == totalSequenceLength);

            uint Bj_2 = H.getTypeOccurrenceB(chosenHaplotype2);
            assert (Bj_2 > 0);

            H.removeType(chosenHaplotype2);

            double piSMCUpdatedHaplotype = csHMM->piSMC(updatedHaplotype,H,c);
            probabilityVector[locus_len_A + i] = Ai*Bj_2*piSMCUpdatedHaplotype/piSMCJointly(chosenHaplotype,chosenHaplotype2,H,c);

            assert (probabilityVector[locus_len_A + i] > 0);
            H.addType(chosenHaplotype2);
        }

        //coalescence (i,*) with (i,*) and coalescence (i,*) with (i,j)
        if (Ai > 1 or Ci > 0) {

            probabilityVector[locus_len_A + vectorHaplotypesB.size()] = Ai*(Ai-1 + Ci)/no_event_pi;
            assert (probabilityVector[locus_len_A + vectorHaplotypesB.size()] > 0);

        } else {

            probabilityVector[locus_len_A + vectorHaplotypesB.size()] = 0;

        }

         //sampling transition move
        sampledIndex = sampleFromVector(probabilityVector);

        if (sampledIndex < locus_len_A){ //have a mutation on (i,*)

            updatedHaplotype = chosenHaplotype;
            updatedHaplotype.type[sampledIndex] = (updatedHaplotype.type[sampledIndex]+1)%2;

            H.addType(updatedHaplotype);

            double Ak = H.getTypeOccurrenceA(updatedHaplotype); //equivalent to Ak+1 when compared to Hk
            event_constant = Ak;
            event_string = "Mutation";

        } else if ((sampledIndex >= locus_len_A) and (sampledIndex < locus_len_A + vectorHaplotypesB.size())){ //coalescence of (i,*) and (*,j)

            chosenHaplotype2 = vectorHaplotypesB[sampledIndex - locus_len_A];
            updatedHaplotype = combineHaplotypes_AandB(chosenHaplotype, chosenHaplotype2); //(i,j)
            assert (updatedHaplotype.type.size() == totalSequenceLength);

            H.removeType(chosenHaplotype2);
            H.addType(updatedHaplotype);

            double Cij = H.getTypeOccurrenceC(updatedHaplotype); //equivalent to Cij+1 when compared to Hk

            assert (Cij > 0);

            // event_constant = 2*n*Cij;
            event_constant = A*B*Cij/(C+1.0);

            event_string = "Coalescence";

        } else if (sampledIndex == locus_len_A+vectorHaplotypesB.size()){ //coalescence of 2 (i,*) or coalescence of (i,*) with (i,j)

            assert (Ai > 0 or Ci > 0);
            
            // event_constant = n*(Ai + 2*Ci - 1);
            event_constant = A*(Ai + Ci - 1);

            event_string = "Coalescence";

        } else {

            assert (1 == 0);

        }

    }

    //if sequence of the form (*,j) is chosen
    if (chosenHaplotype.gamma == "B"){

        assert (Bj > 0);
        vector<Haplotype> vectorHaplotypesA = H.getHaplotypesA();
        probabilityVector.resize(locus_len_B + vectorHaplotypesA.size() + 1, 0);

        //mutation to (*,l)
        for (uint i = 0; i < locus_len_B; i++)
        {
            assert (chosenHaplotype.type.size() == locus_len_B);
            updatedHaplotype = chosenHaplotype;
            updatedHaplotype.type[i] = (updatedHaplotype.type[i]+1)%2;

            double piSMCSampledMutation = csHMM->piSMC(updatedHaplotype, H,c);
            probabilityVector[i] = driving_theta*Bj*piSMCSampledMutation/no_event_pi;
            assert (probabilityVector[i] >= 0);

        }

        //coalescence (*,j) with (i,*)
        for (uint i = 0; i < vectorHaplotypesA.size(); i++){

            chosenHaplotype2 = vectorHaplotypesA[i];
            updatedHaplotype = combineHaplotypes_AandB(chosenHaplotype2, chosenHaplotype); //input (i,*), (*,j) and return (i,j)
            assert (updatedHaplotype.type.size() == totalSequenceLength);

            uint Ai_2 = H.getTypeOccurrenceA(chosenHaplotype2);
            assert (Ai_2 > 0);

            H.removeType(chosenHaplotype2);
            
            double piSMCUpdatedHaplotype = csHMM->piSMC(updatedHaplotype,H,c);
            probabilityVector[locus_len_B + i] = Ai_2*Bj*piSMCUpdatedHaplotype/piSMCJointly(chosenHaplotype,chosenHaplotype2,H,c);
            assert (probabilityVector[locus_len_B + i] > 0);

            H.addType(chosenHaplotype2);
        }

        //coalescence (i,*) with (i,*) and coalescence (i,*) with (i,j)
        if (Bj > 1 or Cj > 0) {

            probabilityVector[locus_len_B + vectorHaplotypesA.size()] = Bj*(Bj-1 + Cj)/no_event_pi;
            assert (probabilityVector[locus_len_B + vectorHaplotypesA.size()] > 0);

        } else {

            probabilityVector[locus_len_B + vectorHaplotypesA.size()] = 0;

        }

         //sampling transition move
        sampledIndex = sampleFromVector(probabilityVector);

        if (sampledIndex < locus_len_B){ //have a mutation on (*,j)
            updatedHaplotype = chosenHaplotype;
            updatedHaplotype.type[sampledIndex] = (updatedHaplotype.type[sampledIndex]+1)%2;

            H.addType(updatedHaplotype);

            double Bl = H.getTypeOccurrenceB(updatedHaplotype); //equivalent to Bl+1 when compared to Hk
            event_constant = Bl;
            event_string = "Mutation";

        } else if ((sampledIndex >= locus_len_B) and (sampledIndex < locus_len_B + vectorHaplotypesA.size())){ //coalescence of (*,j) and (i,*)

            chosenHaplotype2 = vectorHaplotypesA[sampledIndex - locus_len_B];//(i,*)
            updatedHaplotype = combineHaplotypes_AandB(chosenHaplotype2, chosenHaplotype); //(i,j)
            assert (updatedHaplotype.type.size() == totalSequenceLength);

            H.removeType(chosenHaplotype2);
            H.addType(updatedHaplotype);

            double Cij = H.getTypeOccurrenceC(updatedHaplotype); //equivalent to Cij+1 when compared to Hk

            assert (Cij > 0);

            // event_constant = 2*n*Cij;
            event_constant = A*B*Cij/(C+1);

            event_string = "Coalescence";

        } else if (sampledIndex == locus_len_B+vectorHaplotypesA.size()){ //coalescence of 2 (*,j) or coalescence of (*,j) with (i,j)

            assert (Bj > 0 or Cj > 0);

            // event_constant = n*(Bj + 2*Cj - 1);
            event_constant = B*(Bj + Cj - 1);

            event_string = "Coalescence";

        } else {

            assert (1 == 0);

        }
    }

    assert (vectorSum(probabilityVector) > 0);
    assert (probabilityVector[sampledIndex] > 0);
    assert (event_constant > 0);

    double normalisedBackwardTransitionProbability_log = log(probabilityVector[sampledIndex])-log(vectorSum(probabilityVector));

    if (postCollection == false ){

        if (event_string == "Coalescence") {

            for (uint i=0; i < target_mu.size(); i++){

                for (uint j=0; j < target_r.size(); j++){

                    for (uint z=0; z < target_lambda.size(); z++){


                        double N = target_lambda[z]*V;
                        double theta = 4*target_mu[i]*target_lambda[z]*V;
                        double rho = 4*target_r[j]*target_lambda[z]*V;
                        double D = n*(n-1)+(A+C)*locus_len_A*theta+(B+C)*locus_len_B*theta+rho*C;

                        // likelihood_matrix.matrix.at(i).at(j).at(z) += log(event_constant/D);
                        // likelihood_matrix.matrix.at(i).at(j).at(z) += fEventTime(waitingTime, D/(2*N*tau));

                        likelihood_matrix->addValueToElement(log(event_constant/D) + fEventTime(waitingTime, D/(2*N*tau)), i, j, z);
                    }
                }
            }

        } else if (event_string == "Mutation") {

            for (uint i=0; i < target_mu.size(); i++){

                for (uint j=0; j < target_r.size(); j++){

                    for (uint z=0; z < target_lambda.size(); z++){

                        double N = target_lambda[z]*V;
                        double theta = 4*target_mu[i]*target_lambda[z]*V;
                        double rho = 4*target_r[j]*target_lambda[z]*V;
                        double D = n*(n-1)+(A+C)*locus_len_A*theta+(B+C)*locus_len_B*theta+rho*C;
                        
                        // likelihood_matrix.matrix.at(i).at(j).at(z) += log(theta*event_constant/D);
                        // likelihood_matrix.matrix.at(i).at(j).at(z) += fEventTime(waitingTime, D/(2*N*tau));

                        likelihood_matrix->addValueToElement(log(theta*event_constant/D) + fEventTime(waitingTime, D/(2*N*tau)), i, j, z);
                    }
                }
            }

        } else if (event_string == "Recombination") {

            for (uint i=0; i < target_mu.size(); i++){

                for (uint j=0; j < target_r.size(); j++){

                    for (uint z=0; z < target_lambda.size(); z++){

                        double N = target_lambda[z]*V;
                        double theta = 4*target_mu[i]*target_lambda[z]*V;
                        double rho = 4*target_r[j]*target_lambda[z]*V;
                        double D = n*(n-1)+(A+C)*locus_len_A*theta+(B+C)*locus_len_B*theta+rho*C;
                        
                        // likelihood_matrix.matrix.at(i).at(j).at(z) += log(rho*event_constant/D);
                        // likelihood_matrix.matrix.at(i).at(j).at(z) += fEventTime(waitingTime, D/(2*N*tau));

                        likelihood_matrix->addValueToElement(log(rho*event_constant/D) + fEventTime(waitingTime, D/(2*N*tau)), i, j, z);
                    }
                }
            }

        } else {

            assert (1 == 0);
        }

    } else {

        //post collection is true..
        postCollection = false;

        if (event_string == "Coalescence") {

            for (uint i=0; i < target_mu.size(); i++){

                for (uint j=0; j < target_r.size(); j++){

                    for (uint z=0; z < target_lambda.size(); z++){


                        double N = target_lambda[z]*V;
                        double theta = 4*target_mu[i]*target_lambda[z]*V;
                        double rho = 4*target_r[j]*target_lambda[z]*V;
                        double D = n*(n-1)+(A+C)*locus_len_A*theta+(B+C)*locus_len_B*theta+rho*C;
                        
                        // likelihood_matrix.matrix.at(i).at(j).at(z) += log(event_constant/D);
                        // likelihood_matrix.matrix.at(i).at(j).at(z) += probabilityNoEvent(waitingTime, D/(2*N*tau));      

                        likelihood_matrix->addValueToElement(log(event_constant/D) + probabilityNoEvent(waitingTime, D/(2*N*tau)), i, j, z);
                    }
                }
            }

        } else if (event_string == "Mutation") {

            for (uint i=0; i < target_mu.size(); i++){

                for (uint j=0; j < target_r.size(); j++){

                    for (uint z=0; z < target_lambda.size(); z++){

                        double N = target_lambda[z]*V;
                        double theta = 4*target_mu[i]*target_lambda[z]*V;
                        double rho = 4*target_r[j]*target_lambda[z]*V;
                        double D = n*(n-1)+(A+C)*locus_len_A*theta+(B+C)*locus_len_B*theta+rho*C;
                        
                        // likelihood_matrix.matrix.at(i).at(j).at(z) += log(theta*event_constant/D);
                        // likelihood_matrix.matrix.at(i).at(j).at(z) += probabilityNoEvent(waitingTime, D/(2*N*tau));

                        likelihood_matrix->addValueToElement(log(theta*event_constant/D) + probabilityNoEvent(waitingTime, D/(2*N*tau)), i, j, z);
                    }
                }
            }

        } else if (event_string == "Recombination") {

            for (uint i=0; i < target_mu.size(); i++){

                for (uint j=0; j < target_r.size(); j++){

                    for (uint z=0; z < target_lambda.size(); z++){

                        double N = target_lambda[z]*V;
                        double theta = 4*target_mu[i]*target_lambda[z]*V;
                        double rho = 4*target_r[j]*target_lambda[z]*V;
                        double D = n*(n-1)+(A+C)*locus_len_A*theta+(B+C)*locus_len_B*theta+rho*C;

                        // likelihood_matrix.matrix.at(i).at(j).at(z) += log(rho*event_constant/D);
                        // likelihood_matrix.matrix.at(i).at(j).at(z) += probabilityNoEvent(waitingTime, D/(2*N*tau));

                        likelihood_matrix->addValueToElement(log(rho*event_constant/D) + probabilityNoEvent(waitingTime, D/(2*N*tau)), i, j, z);
                    }
                }
            }

        } else {

            assert (1 == 0);
        }
    }

    //just had the last coalescent event - don't need to calculate the corresponding forward probability
    if ((afterLastCollection == true) && (H.getLineageCount() == 1)){
        
        for (uint i=0; i < target_mu.size(); i++){

            for (uint j=0; j < target_r.size(); j++){

                for (uint z=0; z < target_lambda.size(); z++){

                    double N = target_lambda[z]*V;
                    double theta = 4*target_mu[i]*target_lambda[z]*V;
                    double rho = 4*target_r[j]*target_lambda[z]*V;
                    double D = n*(n-1)+(A+C)*locus_len_A*theta+(B+C)*locus_len_B*theta+rho*C;
                    
                    // likelihood_matrix.matrix.at(i).at(j).at(z) -= log(event_constant/D);

                    likelihood_matrix->addValueToElement(-log(event_constant/D), i, j, z);
                }
            }
        }
    }

    assert (normalisedBackwardTransitionProbability_log <= 0);
    return normalisedBackwardTransitionProbability_log;
}

uint ProposalEngine::sampleFromVector(vector<double> vec) {

    double random = sample_uniform_01(*mt_prng) * vectorSum(vec);

    uint counter=0;
    double cumSum = vec[0];

    while (cumSum < random) {

        counter++;
        cumSum += vec[counter];
    }

    assert (counter < vec.size());

    return counter;

}

double ProposalEngine::vectorSum(vector<double> vec) {

    double sum=0;

    for (uint i = 0; i < vec.size(); i++) {

        sum += vec[i];
    }

    return sum;

}

Haplotype ProposalEngine::combineHaplotypes_AandB(Haplotype A, Haplotype B){

    Haplotype C;
    C.gamma = "C";

    vector <bool> C_type(A.type.begin(), A.type.end());
    C_type.insert(C_type.end(), B.type.begin(), B.type.end());
    C.type = C_type;

    return C;

}

Haplotype ProposalEngine::changeHaplotypeCtoA (Haplotype C, uint locusALength){

    vector<bool> A_type(C.type.begin(), C.type.begin() + locusALength);

    Haplotype A;
    A.gamma = "A";
    A.type = A_type;

    return A;

}

Haplotype ProposalEngine::changeHaplotypeCtoB (Haplotype C, uint locusALength){

    vector<bool> B_type(C.type.begin() + locusALength, C.type.end());

    Haplotype B;
    B.gamma = "B";
    B.type = B_type;

    return B;

}

pair <Haplotype,Haplotype> ProposalEngine::splitHaplotypeC (Haplotype C, uint locusALength){

    vector<bool> A_type(C.type.begin(), C.type.begin() + locusALength); //locusA
    vector<bool> B_type(C.type.begin() + locusALength, C.type.end()); //locus B

    Haplotype A;
    Haplotype B;

    A.gamma = "A";
    A.type = A_type;
    B.gamma = "B";
    B.type = B_type;

    return pair<Haplotype,Haplotype> (A,B);

}


uint ProposalEngine::calculateLikelihood(SamplerOutput<LikelihoodSamplesType> * samples_output_local, ModelParameters model_parameters, vector<AncestralTypeMaps> type_data, vector<double> timesCollection, vector<double> viralLoad, pair<uint, uint> seq_lens) {

    double driving_mu = model_parameters.driving_mu;
    double driving_r = model_parameters.driving_rho;
    double driving_lambda = model_parameters.driving_lambda;

    vector<double> target_mu = model_parameters.target_mu;
    vector<double> target_r = model_parameters.target_rho;
    vector<double> target_lambda = model_parameters.target_lambda;

    // Return values stored in class Matrix3d
    Matrix3d<double> likelihood_matrix(target_mu.size(), target_r.size(), target_lambda.size(), 0);

    vector<double> lineages_left(viralLoad.size());

    double V = viralLoad[0];
    double driving_N = driving_lambda*V;
    double driving_theta = 4*driving_mu*driving_N;
    double driving_rho = 4*driving_r*driving_N;

    double driving_theta_A = driving_theta * seq_lens.first;
    double driving_theta_B = driving_theta * seq_lens.second;

    double tau = model_parameters.tau;

    TypeContainer H(mt_prng, type_data[0], seq_lens);
    TypeContainer H_beforeCollection(mt_prng, seq_lens);

    double q_log = 0;

    pair<Haplotype, vector<double> > sampled_haplotype = H.sampleType(driving_theta_A, driving_theta_B, driving_rho);
    Haplotype chosenHaplotype = sampled_haplotype.first;
    double chosenHaplotype_prob = sampled_haplotype.second[0];
    // double backwardEventTimeRate = sampled_haplotype.second[1]/(2*driving_N*tau);

    uint n_t = H.getLineageCount();
    uint C_t = H.getLineageCountC();
    uint A_t = H.getLineageCountA();
    uint B_t = H.getLineageCountB();

    double D = n_t*(n_t-1)+(A_t+C_t)*driving_theta_A+(B_t+C_t)*driving_theta_B+driving_rho*C_t;
    double backwardEventTimeRate = D/(2*driving_N*tau);

    //sample backward event time
    double waitingTime = sampleEventTime(backwardEventTimeRate);
    double waitingTimeSum = waitingTime;
    assert (waitingTime > 0);

    double excessTime = 0;

    bool postCollection = true;
    bool afterLastCollection = false;


    for (uint c = 1; c < timesCollection.size(); c++) {

        while (waitingTimeSum < timesCollection[c]) {

            q_log += fEventTime(waitingTime, backwardEventTimeRate);
            q_log += log(chosenHaplotype_prob);
            q_log += sample_backwardMove_forwardTransition(driving_N, driving_theta, driving_rho, tau, H,
                                            &likelihood_matrix, target_mu, target_r, target_lambda,
                                            chosenHaplotype, waitingTime, c-1, V, postCollection, afterLastCollection);


            // Generate waiting time
            sampled_haplotype = H.sampleType(driving_theta_A, driving_theta_B, driving_rho);
            chosenHaplotype = sampled_haplotype.first;
            chosenHaplotype_prob = sampled_haplotype.second[0];
            // backwardEventTimeRate = sampled_haplotype.second[1]/(2*driving_N*tau);

            uint n_t = H.getLineageCount();
            uint C_t = H.getLineageCountC();
            uint A_t = H.getLineageCountA();
            uint B_t = H.getLineageCountB();

            double D = n_t*(n_t-1)+(A_t+C_t)*driving_theta_A+(B_t+C_t)*driving_theta_B+driving_rho*C_t;
            backwardEventTimeRate = D/(2*driving_N*tau);

            waitingTime = sampleEventTime(backwardEventTimeRate);
            waitingTimeSum += waitingTime;
            assert (waitingTime > 0);

        }

        // have exit the while loop
        excessTime = timesCollection[c] - (waitingTimeSum-waitingTime);
        waitingTimeSum = timesCollection[c];
        assert (excessTime > 0);

        // backward q
        q_log += probabilityNoEvent(excessTime, backwardEventTimeRate);

        // forward transition probabilities
        uint n = H.getLineageCount();
        uint C = H.getLineageCountC();
        uint A = H.getLineageCountA();
        uint B = H.getLineageCountB();

        assert (n == A + B + C);
        
        for (uint i=0; i < target_mu.size(); i++){

            for (uint j=0; j < target_r.size(); j++){

                for (uint z=0; z < target_lambda.size(); z++){

                    double N = target_lambda[z]*V;
                    double theta = 4*target_mu[i]*target_lambda[z]*V;
                    double rho = 4*target_r[j]*target_lambda[z]*V;
                    double D = n*(n-1)+(A+C)*theta*seq_lens.first+(B+C)*theta*seq_lens.second+rho*C;
                    double forwardEventTimeRate = D/(2*N*tau);

                    // likelihood_matrix.matrix.at(i).at(j).at(z) += fEventTime(excessTime, forwardEventTimeRate);

                    likelihood_matrix.addValueToElement(fEventTime(excessTime, forwardEventTimeRate), i, j, z);
                }
            }
        }

        //update H
        H_beforeCollection = H;
        lineages_left[c-1] = n;
        H.addTypeMap(type_data[c]); //ADD NEW AncestralTypeMaps
        postCollection = true;

        // forward transition probabilities: combinatoric terms
        likelihood_matrix.addConstant(probabilityCombinatorics(H, H_beforeCollection));

        //update viral load, driving_N, and driving_theta
        V = viralLoad[c];
        driving_N = driving_lambda*V;
        driving_theta = 4*driving_mu*driving_N;
        driving_rho = 4*driving_r*driving_N;

        driving_theta_A = driving_theta * seq_lens.first;
        driving_theta_B = driving_theta * seq_lens.second;

        // Generate waiting time
        sampled_haplotype = H.sampleType(driving_theta_A, driving_theta_B, driving_rho);
        chosenHaplotype = sampled_haplotype.first;
        chosenHaplotype_prob = sampled_haplotype.second[0];
        // backwardEventTimeRate = sampled_haplotype.second[1]/(2*driving_N*tau);

        uint n_t = H.getLineageCount();
        uint C_t = H.getLineageCountC();
        uint A_t = H.getLineageCountA();
        uint B_t = H.getLineageCountB();

        double D = n_t*(n_t-1)+(A_t+C_t)*driving_theta_A+(B_t+C_t)*driving_theta_B+driving_rho*C_t;
        backwardEventTimeRate = D/(2*driving_N*tau);

        waitingTime = sampleEventTime(backwardEventTimeRate);
        waitingTimeSum += waitingTime;
        assert (waitingTime > 0);
    }

    uint secret_counter = 0;
    uint secret_threshold = 10000;

    //after last collection time:
    afterLastCollection = true;
    n_t = H.getLineageCount();
    C_t = H.getLineageCountC();
    A_t = H.getLineageCountA();
    B_t = H.getLineageCountB();

    while (n_t > 5) {
        q_log += fEventTime(waitingTime, backwardEventTimeRate);
        q_log += log(chosenHaplotype_prob);
        q_log += sample_backwardMove_forwardTransition(driving_N, driving_theta, driving_rho, tau, H,
                                        &likelihood_matrix, target_mu, target_r, target_lambda,
                                        chosenHaplotype, waitingTime, timesCollection.size()-1,V, postCollection, afterLastCollection);

        // Generate waiting time
        sampled_haplotype = H.sampleType(driving_theta_A, driving_theta_B, driving_rho);
        chosenHaplotype = sampled_haplotype.first;
        chosenHaplotype_prob = sampled_haplotype.second[0];
        // backwardEventTimeRate = sampled_haplotype.second[1]/(2*driving_N*tau);

        n_t = H.getLineageCount();
        C_t = H.getLineageCountC();
        A_t = H.getLineageCountA();
        B_t = H.getLineageCountB();

        double D = n_t*(n_t-1)+(A_t+C_t)*driving_theta_A+(B_t+C_t)*driving_theta_B+driving_rho*C_t;
        backwardEventTimeRate = D/(2*driving_N*tau);

        waitingTime = sampleEventTime(backwardEventTimeRate);
        assert (waitingTime > 0);
        waitingTimeSum += waitingTime;
    }

    lineages_left[timesCollection.size()-1] = n_t;
    waitingTimeSum -= waitingTime;

    assert(waitingTimeSum > 0);

    double permutation_constant = 0;

    //add permutation terms
    for (int i=0; i < timesCollection.size(); i++) {

        permutation_constant += probabilitySamplePermutation(type_data[i]);
    }

    assert(permutation_constant < 0);
    //assert(q_log < 0);

    Matrix3d<double> likelihood_matrix_before_final_update = likelihood_matrix;

    // subtract forward probabilities by backward probability and add probability of observing the data (permutation_constant)
    likelihood_matrix.addConstant(permutation_constant - q_log);
	
	// subtract forward probabilities by the uniform distributions for the remaining lineages
	likelihood_matrix.addConstant( - (double)(seq_lens.first + seq_lens.second)*log(2)*C_t - (double)(seq_lens.first)*log(2)*A_t - (double)(seq_lens.second)*log(2)*B_t);

    for (uint i = 0; i < likelihood_matrix.dim1(); i++) {

        for (uint j = 0; j < likelihood_matrix.dim2(); j++) {

            for (uint k = 0; k < likelihood_matrix.dim3(); k++) {

                if (likelihood_matrix.getElement(i,j,k) >= 0) {

                    cout << "WARNING: non negative combined log likelihood encountered " << likelihood_matrix.getElement(i,j,k) << endl; 
                    cout << "\t for i=" << i << ", j=" << j << ", k=" << k << endl;
                    cout << "\t for mu=" << target_mu.at(i) << ", rho=" << target_r.at(j) << ", lambda=" << target_lambda.at(k) << endl;
                    cout << "\t likelihood_matrix_before_final_update " << likelihood_matrix_before_final_update.getElement(i,j,k) << endl;
                    cout << "\t permutation_constant " << permutation_constant << endl;
                    cout << "\t P(mrca) term " << (double)(seq_lens.first + seq_lens.second)*log(2) << endl;
                    cout << "\t q_log" << q_log << endl;
					cout << "\t waitingTimeSum" << waitingTimeSum << endl;
					cout << "\t for C_t=" << C_t << ", A_t=" << A_t << ", B_t=" << B_t << endl;
					
                }
            }
        }
    }

    samples_output_local->addSample(likelihood_matrix, waitingTimeSum, lineages_left);

    return 1;
}



