
/*
ImportanceSampler.cpp - This file is part of the Coalescenator (v1.0.0)


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


#include <math.h>
#include <sstream>
#include <random>
#include <chrono>

#include "ImportanceSampler.hpp"
#include "Utils.hpp"
#include "Matrix3d.hpp"
#include "ConditionalSamplingHMM.hpp"
#include "ProposalEngine.hpp"
#include "LikelihoodSamples.hpp"


void ImportanceSampler::samplerCallback(Sample * sample, SamplerOutput<LikelihoodSamplesType> * sampler_output_global, ModelParameters model_parameters, uint thread_iterations, uint quadrature_points, uint seed, bool deterministic) {

	unique_lock<mutex> global_lock(sampler_mutex);

	const vector<AncestralTypeMaps> type_data = sample->getTypeMaps();
	const vector<double> collection_times = sample->getTimes(); 
	const vector<double> viral_loads = sample->getViralLoads(); 
	const pair <uint, uint> seq_lens = sample->getSeqLens(); 

	global_lock.unlock();
	
	assert (type_data.size() == viral_loads.size());
	assert (collection_times.size() == viral_loads.size());

	auto thread_id = this_thread::get_id();

	// Seen PRNG using thead ids to ensure somewhat independent seeds for each thread
	uint prng_seed = seed;

	if (!deterministic) {

		prng_seed *= hash<thread::id>()(thread_id);
	}

	mt19937 mt_prng(prng_seed);

	vector <double> scaled_driving_theta(collection_times.size(), 0);
	vector <double> scaled_driving_rho(collection_times.size(), 0);

	for (uint i = 0; i < collection_times.size(); i++) {

		scaled_driving_theta[i] = 4 * model_parameters.driving_lambda * viral_loads[i] * model_parameters.driving_mu;
    	scaled_driving_rho[i] = 4 * model_parameters.driving_lambda * viral_loads[i] * model_parameters.driving_rho;
	}

	uint max_seq_len;

	if (seq_lens.first > seq_lens.second) {

		max_seq_len = seq_lens.first;

	} else {

		max_seq_len = seq_lens.second;
	}

	SamplerOutput<LikelihoodSamplesType> sampler_output_local(model_parameters.target_mu.size(), model_parameters.target_rho.size(), model_parameters.target_lambda.size(), thread_iterations, collection_times.size());
	ConditionalSamplingHMM csHMM(quadrature_points, scaled_driving_rho, scaled_driving_theta, max_seq_len);
	ProposalEngine proposal_engine(&mt_prng, &csHMM);

	uint sample_counter = 0;

	while (sample_counter < thread_iterations) {
		
				// Update local matrix
		sample_counter += proposal_engine.calculateLikelihood(&sampler_output_local, model_parameters, type_data, collection_times, viral_loads, seq_lens);
	}
	
	assert (sample_counter == thread_iterations);	

	// Update global matrices
	global_lock.lock();
	sampler_output_global->addSamples(sampler_output_local);
	global_lock.unlock();
}


SamplerOutput<double> ImportanceSampler::runSampler(Sample * sample, ModelParameters model_parameters, uint num_iterations, uint num_threads, uint quadrature_points, uint seed) {

	assert (num_threads > 0);
	uint thread_iterations = floor(num_iterations/num_threads);
	uint remaining_iterations = num_iterations % num_threads;

	SamplerOutput<LikelihoodSamplesType> sampler_output(model_parameters.target_mu.size(), model_parameters.target_rho.size(), model_parameters.target_lambda.size(), num_iterations, sample->getTimes().size());

	bool deterministic = false;

	if (num_threads == 1) {

		deterministic = true;
	}

	vector<thread> threads; 

	// Launch first thread
	threads.push_back(thread(&ImportanceSampler::samplerCallback, this, sample, &sampler_output, model_parameters, thread_iterations + remaining_iterations, quadrature_points, seed, deterministic));

	// Launch threads
	for (uint i = 0; i < (num_threads-1); i++) {
	
		threads.push_back(thread(&ImportanceSampler::samplerCallback, this, sample, &sampler_output, model_parameters, thread_iterations, quadrature_points, seed, deterministic));
	}	

   for(auto& thread : threads) {

        thread.join();
    }

    SamplerOutput<double> results(model_parameters.target_mu.size(), model_parameters.target_rho.size(), model_parameters.target_lambda.size(), 0, sample->getTimes().size());
	results.sumSamples(sampler_output);

	return results;
}


