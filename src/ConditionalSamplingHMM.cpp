
/*
ConditionalSamplingHMM.cpp - This file is part of the Coalescenator (v1.0.0)


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


#include "ConditionalSamplingHMM.hpp"

#include <math.h>

#include "boost/math/special_functions/binomial.hpp"

ConditionalSamplingHMM::ConditionalSamplingHMM(uint n_quadrature_points, vector<double> & r_rates_in, vector<double> & thetas_in, uint max_seq_len) {

	// Laguerre-Gauss quadrature points from Table 25.9, Handbook of Mathematical Functions, Abramowitz and Stegun, 1972
	if (n_quadrature_points == 4) {

		double quadrature4[] = {0, 0.415774556783, 2.294280360279, 6.289945082937, 1E+300};
		quadratures = vector<double> (quadrature4, quadrature4 + sizeof(quadrature4)/sizeof(double));

	} else if (n_quadrature_points == 8) {

		double quadrature8[] = {0, 0.193043676560, 1.026664895339, 2.567876744951, 4.900353084526, 8.182153444563, 12.7344180291798, 19.3957278622, 1E+300};
		quadratures = vector<double> (quadrature8, quadrature8 + sizeof(quadrature8)/sizeof(double));

	} else if (n_quadrature_points == 16) {

		double quadrature16[] = {0, 0.093307812017, 0.492691740302, 1.215595412071, 2.269949526204, 3.667622721751, 5.425336627414, 7.565916226613, 10.120228568019, 13.130282482176, 16.6544077083, 20.7764788994, 25.6238942267, 31.407519169754, 38.530683306486, 48.026085572686, 1E+300};
		quadratures = vector<double> (quadrature16, quadrature16 + sizeof(quadrature16)/sizeof(double));

	} else {

		cerr << "ERROR: Request for " << n_quadrature_points << " not supported. Only n_quadrature_points in {4,8,16} supported." << endl;
		exit(1);
	}

	r_rates = r_rates_in;
	thetas = thetas_in;

	vector <double> weights_temp(quadratures.size() - 1);

	for (uint i = 1; i < quadratures.size(); i++) {

		weights_temp[i - 1] = exp(-quadratures[i-1]) - exp(-quadratures[i]);
		assert(weights_temp[i - 1] > 0);
	}

	weights = weights_temp;

	e_binomial = Eigen::ArrayXXd::Zero(max_seq_len + 1, max_seq_len + 1);

	for (uint i = 0; i < max_seq_len + 1; i++) {

		for (uint j = i; j < max_seq_len + 1; j++) {

			e_binomial(j,i) = boost::math::binomial_coefficient<double>(j,i);
			assert (e_binomial(j,i) > 0);
		}
	}
}

void ConditionalSamplingHMM::calcConstantsTransition(uint n_all, uint time_idx) {

	assert(r_rates.size() > time_idx);

	HMMConstants_transition hmm_constants_transition;
	hmm_constants_transition.y = vector <double>(weights.size());
	hmm_constants_transition.z = Eigen::ArrayXXd::Zero(weights.size(), weights.size());

	double n_r_constant1 = n_all / (r_rates[time_idx] + n_all);

	if (double_compare(r_rates[time_idx], n_all)) {

		for (uint i = 0; i < weights.size(); i++) {

			hmm_constants_transition.y[i] = n_r_constant1 / weights[i] * (exp(-quadratures[i]/n_r_constant1) - exp(-quadratures[i+1]/n_r_constant1));

			if (hmm_constants_transition.y[i] > 1) {

				cerr << "y in transition probability is larger than one: " << hmm_constants_transition.y[i] << endl;
				assert(false);
			}

			assert (hmm_constants_transition.y[i] >= 0);
			assert (hmm_constants_transition.y[i] <= 1);

			for (uint j = 0; j < weights.size(); j++) {

				if (i > j) {

					hmm_constants_transition.z(i,j) =  (weights[j] + (quadratures[j] * exp(-quadratures[j]) - quadratures[j+1] * exp(-quadratures[j+1])));

				} else if (i < j) {

					hmm_constants_transition.z(i,j) = weights[j] / weights[i] * (weights[i] + (quadratures[i] * exp(-quadratures[i]) - quadratures[i+1] * exp(-quadratures[i+1])));

				} else {

					double term1 = weights[i] * (weights[i] + (quadratures[i] * exp(-quadratures[i]) - quadratures[i+1] * exp(-quadratures[i+1])));
					double term2 =  - (quadratures[i] - quadratures[i+1]) * exp(-(quadratures[i] + quadratures[i+1])) - 1/2 * (exp(-2*quadratures[i]) - exp(-2*quadratures[i+1]));

					hmm_constants_transition.z(i,j) = 1 / weights[i] * (term1 + term2);
				}

				assert (hmm_constants_transition.z(i,j) >= 0);
				assert (hmm_constants_transition.z(i,j) <= 1);
			}
		}

	} else {

		double n_r_constant2 = r_rates[time_idx] / (r_rates[time_idx] - n_all);
		double n_r_constant3 = n_all / r_rates[time_idx];

		for (uint i = 0; i < weights.size(); i++) {

			hmm_constants_transition.y[i] = n_r_constant1 / weights[i] * (exp(-quadratures[i]/n_r_constant1) - exp(-quadratures[i+1]/n_r_constant1));

			if (hmm_constants_transition.y[i] > 1) {

				cerr << "y in transition probability is larger than one: " << hmm_constants_transition.y[i] << endl;
				assert(false);
			}


			assert (hmm_constants_transition.y[i] >= 0);
			assert (hmm_constants_transition.y[i] <= 1);

			for (uint j = 0; j < weights.size(); j++) {

				if (i > j) {

					hmm_constants_transition.z(i,j) = n_r_constant2 * (weights[j] - n_r_constant3 * (exp(-quadratures[j]/n_r_constant3) - exp(-quadratures[j+1]/n_r_constant3)));

				} else if (i < j) {

					hmm_constants_transition.z(i,j) = n_r_constant2 * weights[j] / weights[i] * (weights[i] - n_r_constant3 * (exp(-quadratures[i]/n_r_constant3) - exp(-quadratures[i+1]/n_r_constant3)));

				} else {

					double term1 = weights[i] * (weights[i] - n_r_constant3 * (exp(-quadratures[i]/n_r_constant3) - exp(-quadratures[i+1]/n_r_constant3)));
					double term2 = n_r_constant1 / n_r_constant2 * (exp(-quadratures[i]/n_r_constant1) - exp(-quadratures[i+1]/n_r_constant1));
					double term3 = n_r_constant3 * (exp(-quadratures[i])*exp(-quadratures[i+1]/n_r_constant3) - exp(-quadratures[i+1])*exp(-quadratures[i]/n_r_constant3));

					hmm_constants_transition.z(i,j) = n_r_constant2 / weights[i] * (term1 - term2 - term3);
				}

				assert (hmm_constants_transition.z(i,j) >= 0);
				assert (hmm_constants_transition.z(i,j) <= 1);
			}
		}
	}

	vector <uint> key_vec(2);
	key_vec[0] = n_all;
	key_vec[1] = time_idx;

	constant_container_transition[key_vec] = hmm_constants_transition;
}

void ConditionalSamplingHMM::calcConstantsEmission(uint n_all, uint time_idx, uint seq_len) {

	Eigen::ArrayXXd inner_map = Eigen::ArrayXXd::Zero(weights.size(), (seq_len + 1));

	assert(thetas.size() > time_idx);

	for (uint i = 1; i < quadratures.size(); i++) {

		for (uint j = 0; j <= seq_len; j++) {  // Loop over number of differences between query and absorbing haplotype

    		double e_prob = 0;

	    	for (uint l = 0; l <= j; l++) {

				double inner_sum = 0;

				for (uint k = 0; k <= (seq_len - j); k++) {

		    		double rate = (2*thetas[time_idx]*(k+l)/n_all) +  1;
		    		assert (rate > 0);

		    		inner_sum += e_binomial(seq_len - j, k) / rate * (exp(-rate * quadratures[i - 1]) - exp(-rate * quadratures[i]));
		       	}

		    	e_prob += inner_sum * e_binomial(j, l) * pow(-1,(double)l);
		    }

			e_prob /= (weights[i-1] * pow(2,(double) seq_len));
            e_prob = max(e_prob, double_precision);
			inner_map(i-1,j) = e_prob;

			assert (e_prob >= 0);
			assert (e_prob <= 1);
		}
	}

	vector <uint> key_vec_outer(3);
	key_vec_outer[0] = n_all;
	key_vec_outer[1] = time_idx;
	key_vec_outer[2] = seq_len;

	constant_container_emission[key_vec_outer] = inner_map;
}

double ConditionalSamplingHMM::initialDensity(uint n_all, uint n_type, double weight) {

	return n_type / (double) n_all * weight;
}

// Simple handling of emission from unobserved
double ConditionalSamplingHMM::calculatePartialHaplotypeEmissionProbability(uint discretisation_idx, vector<uint> & key_vec_emission, vector<int> & haplotype_counts, vector<int> & difference_vector) {

	uint num_lineages_complete = 0;
	double impute_prob = 0;

	for (uint i = 0; i < haplotype_counts.size(); i++) {

		if (difference_vector[i] >= 0) { // If complete haplotype

			impute_prob += haplotype_counts[i] * constant_container_emission[key_vec_emission](discretisation_idx, difference_vector[i]);
			num_lineages_complete += haplotype_counts[i];
			assert (haplotype_counts[i] > 0);
		}
	}

	if (num_lineages_complete == 0) {

		impute_prob = 1/ pow(2,(double) key_vec_emission[2]);

	} else {

		impute_prob /= num_lineages_complete;

	}

	assert(impute_prob >= 0);
	assert(impute_prob <= 1);

	return impute_prob;
}

double ConditionalSamplingHMM::piSMC(Haplotype & cur_type, TypeContainer & type_container, uint time_idx) {

	uint cur_lineages = type_container.getLineageCount();
	uint seq_len_A = type_container.getLocusLengthA();
	uint seq_len_B = type_container.getLocusLengthB();

	vector<vector<int> > differences = type_container.getDifferences(cur_type);

	vector <int> haplotype_counts = differences[0];
	vector <int> difference_A = differences[1];
	vector <int> difference_B = differences[2];

	if (haplotype_counts.size() == 0) {

		assert (difference_A.size() == 0);
		assert (difference_B.size() == 0);
		assert (cur_lineages == 0);

		if (cur_type.gamma == "A") {

			return 1/ pow(2,(double) (seq_len_A));

		} else if (cur_type.gamma == "B") {

			return 1/ pow(2,(double) (seq_len_B));

		} else {

			return 1/ pow(2,(double) (seq_len_A + seq_len_B));

		}


	}

	vector <uint> key_vec_transition(2);
	key_vec_transition[0] = cur_lineages;
	key_vec_transition[1] = time_idx;

	if (constant_container_transition.count(key_vec_transition) == 0) {

		calcConstantsTransition(cur_lineages, time_idx);
	}

	vector <uint> key_vec_emission_A(3);
	key_vec_emission_A[0] = cur_lineages;
	key_vec_emission_A[1] = time_idx;
	key_vec_emission_A[2] = seq_len_A;

	if (constant_container_emission.count(key_vec_emission_A) == 0) {

		calcConstantsEmission(cur_lineages, time_idx, seq_len_A);
	}

	vector <uint> key_vec_emission_B(3);
	key_vec_emission_B[0] = cur_lineages;
	key_vec_emission_B[1] = time_idx;
	key_vec_emission_B[2] = seq_len_B;

	if (constant_container_emission.count(key_vec_emission_B) == 0) {

		calcConstantsEmission(cur_lineages, time_idx, seq_len_B);
	}

	double pi_smc = 0;

	// Partial haplotype probabilities calculated using marginal locus distribution
	if (cur_type.gamma == "A") {

		for (uint i = 0; i < weights.size(); i++) {

			for (uint idx_A = 0; idx_A < difference_A.size(); idx_A++) {

				double initial_density = initialDensity(cur_lineages, haplotype_counts[idx_A], weights[i]);
				assert (initial_density > 0);

				if (difference_A[idx_A] == -1) {

					 pi_smc += initial_density * calculatePartialHaplotypeEmissionProbability(i, key_vec_emission_A, haplotype_counts, difference_A);

				} else {

					 pi_smc += initial_density * constant_container_emission[key_vec_emission_A](i, difference_A[idx_A]);
				}
			}
		}

	} else if (cur_type.gamma == "B") {

		for (uint idx_B = 0; idx_B < difference_B.size(); idx_B++) {

			for (uint j = 0; j < weights.size(); j++) {

				double initial_density = initialDensity(cur_lineages, haplotype_counts[idx_B], weights[j]);
				assert (initial_density > 0);

				if (difference_B[idx_B] == -1) {

					 pi_smc += initial_density * calculatePartialHaplotypeEmissionProbability(j, key_vec_emission_B, haplotype_counts, difference_B);

				} else {

					 pi_smc += initial_density * constant_container_emission[key_vec_emission_B](j, difference_B[idx_B]);
				}
			}
		}

	} else {

		// Init base matrix
		Eigen::ArrayXXd base_matrix = Eigen::ArrayXXd::Zero(weights.size(), difference_A.size());
		Eigen::ArrayXd base_sum_vector = Eigen::ArrayXd::Zero(weights.size());

		for (uint i = 0; i < weights.size(); i++) {

			for (uint idx_A = 0; idx_A < difference_A.size(); idx_A++) {

				// If partial haplotype emission
				if (difference_A[idx_A] == -1) {

					base_matrix(i,idx_A) = calculatePartialHaplotypeEmissionProbability(i, key_vec_emission_A, haplotype_counts, difference_A) * initialDensity(cur_lineages, haplotype_counts[idx_A], weights[i]);
					base_sum_vector(i) += base_matrix(i,idx_A);

				} else {

					base_matrix(i,idx_A) = constant_container_emission[key_vec_emission_A](i, difference_A[idx_A]) * initialDensity(cur_lineages, haplotype_counts[idx_A], weights[i]);
					base_sum_vector(i) += base_matrix(i,idx_A);
				}

				assert (base_matrix(i,idx_A) >= 0);
			}
		}

		for (uint idx_B = 0; idx_B < difference_B.size(); idx_B++) {

			for (uint j = 0; j < weights.size(); j++) {

				double F_2 = 0;

				for (uint i = 0; i < weights.size(); i++) {

					F_2 += base_sum_vector(i) * constant_container_transition[key_vec_transition].z(i,j);
				}

				F_2 *= haplotype_counts[idx_B] / (double) cur_lineages;
				F_2 += constant_container_transition[key_vec_transition].y[j] * base_matrix(j,idx_B);

				// If partial haplotype emission
				if (difference_B[idx_B] == -1) {

					pi_smc += F_2 * calculatePartialHaplotypeEmissionProbability(j, key_vec_emission_B, haplotype_counts, difference_B);

				} else {

					pi_smc += F_2 * constant_container_emission[key_vec_emission_B](j, difference_B[idx_B]);
				}

			}
		}
	}

	// if (double_underflow > pi_smc) {

	// 	if (cur_type.gamma == "A") {

	// 		return log(1/ pow(2,(double) (seq_len_A)));

	// 	} else if (cur_type.gamma == "B") {

	// 		return log(1/ pow(2,(double) (seq_len_B)));

	// 	} else {

	// 		return log(1/ pow(2,(double) (seq_len_A + seq_len_B)));

	// 	}

	// }

	assert (pi_smc >= 0);
	assert (pi_smc <= 1);

	return (pi_smc);
}
