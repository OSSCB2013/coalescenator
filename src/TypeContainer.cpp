
/*
TypeContainer.cpp - This file is part of the Coalescenator (v1.0.0)


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


#include "TypeContainer.hpp"

#include <vector>
#include <unordered_map>
#include <random>
#include <map>

#include "Utils.hpp"


using namespace std;

TypeContainer::TypeContainer (mt19937* mt_rng_pt_in, AncestralTypeMaps & type_map_in, pair<uint, uint> & loci_lengths_in) {

	mt_rng_pt = mt_rng_pt_in;
	uniform_01_sampler_p = new uniform_01_sampler_t(0.0,1.0);

	type_map_C = type_map_in.C;
	type_map_A = type_map_in.A;
	type_map_B = type_map_in.B;

	num_lineages_A = 0;
	num_lineages_B = 0;
	num_lineages_C = 0;

    second_locus_start_idx = loci_lengths_in.first;
    locus_length_A = loci_lengths_in.first;
    locus_length_B = loci_lengths_in.second;
   	haplotype_length = locus_length_A + locus_length_B;

	fetched_overlap = false;

	// Create loci indices
	for (auto it = type_map_C.begin(); it != type_map_C.end(); it++) {

		num_lineages_C += it->second;

		assert (it->second > 0);
		assert (it->first.size() == haplotype_length);

		pair<vector<bool>, vector<bool> > loci_types;
		getLociTypes(it->first, &loci_types);

		assert (it->first.size() == loci_types.first.size() + loci_types.second.size());

		if (!(type_map_Ca.insert({loci_types.first, it->second}).second)) {

			type_map_Ca.at(loci_types.first) += it->second;
		}

		if (!(type_map_Cb.insert({loci_types.second, it->second}).second)) {

			type_map_Cb.at(loci_types.second) += it->second;
		}
	}

	for (auto it = type_map_A.begin(); it != type_map_A.end(); it++) {

		assert (it->second > 0);
		assert (it->first.size() == locus_length_A);
		num_lineages_A += it->second;
	}

	for (auto it = type_map_B.begin(); it != type_map_B.end(); it++) {

		assert (it->second > 0);
		assert (it->first.size() == locus_length_B);
		num_lineages_B += it->second;
	}
}


TypeContainer::TypeContainer(mt19937* mt_rng_pt_in, pair<uint, uint> & loci_lengths_in) {

	mt_rng_pt = mt_rng_pt_in;
	uniform_01_sampler_p = new uniform_01_sampler_t(0.0,1.0);

	num_lineages_A = 0;
	num_lineages_B = 0;
	num_lineages_C = 0;

    second_locus_start_idx = loci_lengths_in.first;
    locus_length_A = loci_lengths_in.first;
    locus_length_B = loci_lengths_in.second;
   	haplotype_length = locus_length_A + locus_length_B;

	fetched_overlap = false;
}


TypeContainer::TypeContainer(const TypeContainer & in_copy) {

    mt_rng_pt = in_copy.mt_rng_pt;
    uniform_01_sampler_p = new uniform_01_sampler_t(0.0,1.0);

	type_map_A = in_copy.type_map_A;
    type_map_B = in_copy.type_map_B;
    type_map_C = in_copy.type_map_C;
    type_map_Ca = in_copy.type_map_Ca;
    type_map_Cb = in_copy.type_map_Cb;

    num_lineages_A = in_copy.num_lineages_A;
    num_lineages_B = in_copy.num_lineages_B;
    num_lineages_C = in_copy.num_lineages_C;

    second_locus_start_idx = in_copy.second_locus_start_idx;
    haplotype_length = in_copy.haplotype_length;
    locus_length_A = in_copy.locus_length_A;
    locus_length_B = in_copy.locus_length_B;

    fetched_overlap = in_copy.fetched_overlap;
    type_map_overlap = in_copy.type_map_overlap;
}

// Non exception safe
TypeContainer & TypeContainer::operator=(const TypeContainer & in_copy) {

	if (this != &in_copy) {

		delete uniform_01_sampler_p;

	    mt_rng_pt = in_copy.mt_rng_pt;
	    uniform_01_sampler_p = new uniform_01_sampler_t(0.0,1.0);

		type_map_A = in_copy.type_map_A;
	    type_map_B = in_copy.type_map_B;
	    type_map_C = in_copy.type_map_C;
	    type_map_Ca = in_copy.type_map_Ca;
	    type_map_Cb = in_copy.type_map_Cb;

	    num_lineages_A = in_copy.num_lineages_A;
	    num_lineages_B = in_copy.num_lineages_B;
	    num_lineages_C = in_copy.num_lineages_C;

	    second_locus_start_idx = in_copy.second_locus_start_idx;
	    haplotype_length = in_copy.haplotype_length;
	    locus_length_A = in_copy.locus_length_A;
	    locus_length_B = in_copy.locus_length_B;

	    fetched_overlap = in_copy.fetched_overlap;
	    type_map_overlap = in_copy.type_map_overlap;
	}

	return *this; 
}


TypeContainer::~TypeContainer() {

	delete uniform_01_sampler_p;
}


void TypeContainer::getLociTypes(const vector<bool> & target_vector, pair<vector<bool>, vector<bool> > * loci_types) {

	loci_types->first.assign(target_vector.begin(), target_vector.begin() + second_locus_start_idx);
	loci_types->second.assign(target_vector.begin() + second_locus_start_idx, target_vector.end());
}


bool TypeContainer::isEmpty() {

	if (num_lineages_A == 0 and num_lineages_B == 0 and num_lineages_C == 0) {

		assert(type_map_A.empty());
		assert(type_map_B.empty());
		assert(type_map_C.empty());
		assert(type_map_Ca.empty());
		assert(type_map_Cb.empty());

		return true;

	} else {

		return false;
	}
}


uint TypeContainer::getLineageCount() {

	uint num_lineages = num_lineages_A + num_lineages_B + num_lineages_C;

	return num_lineages;
}


uint TypeContainer::getLineageCountA() {

	return num_lineages_A;
}


uint TypeContainer::getLineageCountB() {

	return num_lineages_B;
}


uint TypeContainer::getLineageCountC() {

	return num_lineages_C;
}


uint TypeContainer::getLocusLengthA() {

	return locus_length_A;
}


uint TypeContainer::getLocusLengthB() {

	return locus_length_B;
}


uint TypeContainer::getTotalSequenceLength() {

	return haplotype_length;
}


uint TypeContainer::getUniqueCount() {

	uint unique_counts = type_map_A.size() + type_map_B.size() + type_map_C.size();

	return unique_counts;
}


uint TypeContainer::getTypeOccurrenceC(Haplotype & haplotype) {

    if (haplotype.gamma == "C") {

    	if (type_map_C.count(haplotype.type) > 0) {

    		return type_map_C.at(haplotype.type);
    	}
	}

    return 0;
}


uint TypeContainer::getTypeOccurrenceA(Haplotype & haplotype) {

 	if (haplotype.gamma == "C") {

 		pair<vector<bool>, vector<bool> > loci_types;
		getLociTypes(haplotype.type, &loci_types);

 		if (type_map_A.count(loci_types.first) > 0) {

        	return type_map_A.at(loci_types.first);

        } else {

        	return 0;
        }

    } else if (haplotype.gamma == "A") {

    	if (type_map_A.count(haplotype.type) > 0) {

        	return type_map_A.at(haplotype.type);

    	} else {

    		return 0;
    	}
    
    } else {

	    assert(haplotype.gamma == "B");
	    return 0;
	}
}


uint TypeContainer::getTypeOccurrenceB(Haplotype & haplotype) {

 	if (haplotype.gamma == "C"){

 		pair<vector<bool>, vector<bool> > loci_types;
		getLociTypes(haplotype.type, &loci_types);

 		if (type_map_B.count(loci_types.second) > 0) {

        	return type_map_B.at(loci_types.second);

        } else {

        	return 0;
        }

    } else if (haplotype.gamma == "B") {

    	if (type_map_B.count(haplotype.type) > 0) {

        	return type_map_B.at(haplotype.type);

    	} else {

    		return 0;
    	}
    }

	assert(haplotype.gamma == "A");
    return 0;
}


uint TypeContainer::getTypeOccurrenceCa(Haplotype & haplotype) {

    if (haplotype.gamma == "C") {

 		pair<vector<bool>, vector<bool> > loci_types;
		getLociTypes(haplotype.type, &loci_types);

    	if (type_map_Ca.count(loci_types.first) > 0) {

        	return type_map_Ca.at(loci_types.first);

        } else {

        	return 0;
        }

    } else if (haplotype.gamma == "A"){

        if (type_map_Ca.count(haplotype.type) > 0) {

        	return type_map_Ca.at(haplotype.type);

    	} else {

    		return 0;
    	}
    }

    assert(haplotype.gamma == "B");
    return 0;
}


uint TypeContainer::getTypeOccurrenceCb(Haplotype & haplotype) {

  	if (haplotype.gamma == "C") {

 		pair<vector<bool>, vector<bool> > loci_types;
		getLociTypes(haplotype.type, &loci_types);

  		if (type_map_Cb.count(loci_types.second) > 0) {

        	return type_map_Cb.at(loci_types.second);

        } else {

        	return 0;
        }

    } else if (haplotype.gamma == "B") {

        if (type_map_Cb.count(haplotype.type) > 0) {

        	return type_map_Cb.at(haplotype.type);

    	} else {

    		return 0;
    	}
    }

    assert(haplotype.gamma == "A");
    return 0;
}


void TypeContainer::addType(Haplotype & haplotype) {

	if (haplotype.gamma == "A") {

		num_lineages_A++;
		assert (haplotype.type.size() == locus_length_A);

		if (!(type_map_A.insert({haplotype.type, 1}).second)) {

			type_map_A.at(haplotype.type)++;
		} 

	} else if (haplotype.gamma == "B") {

		num_lineages_B++;
		assert (haplotype.type.size() == locus_length_B);

		if (!(type_map_B.insert({haplotype.type, 1}).second)) {

			type_map_B.at(haplotype.type)++;
		} 

	} else {

		num_lineages_C++;
		assert (haplotype.type.size() == haplotype_length);
		assert(haplotype.gamma == "C");

		if (!(type_map_C.insert({haplotype.type, 1}).second)) {

			type_map_C.at(haplotype.type)++;
		} 

		// Update locus counts
 		pair<vector<bool>, vector<bool> > loci_types;
		getLociTypes(haplotype.type, &loci_types);

		if (!(type_map_Ca.insert({loci_types.first, 1}).second)) {

			type_map_Ca.at(loci_types.first)++;
		} 

		if (!(type_map_Cb.insert({loci_types.second, 1}).second)) {

			type_map_Cb.at(loci_types.second)++;
		} 
	}
}


void TypeContainer::removeType(Haplotype & haplotype) {

	if (haplotype.gamma == "A") {

		num_lineages_A--;

		if (type_map_A.at(haplotype.type) > 1) {

			type_map_A.at(haplotype.type)--;

		} else {

			assert(type_map_A.erase(haplotype.type));
		}

	} else if (haplotype.gamma == "B") {

		num_lineages_B--;

		if (type_map_B.at(haplotype.type) > 1) {

			type_map_B.at(haplotype.type)--;

		} else {

			assert(type_map_B.erase(haplotype.type));
		}

	} else {

		num_lineages_C--;

		if (type_map_C.at(haplotype.type) > 1) {

			type_map_C.at(haplotype.type)--;

		} else {

			assert(type_map_C.erase(haplotype.type));
		}

		// Update locus counts
 		pair<vector<bool>, vector<bool> > loci_types;
		getLociTypes(haplotype.type, &loci_types);

		if (type_map_Ca.at(loci_types.first) > 1) {

			type_map_Ca.at(loci_types.first)--;

		} else {

			assert(type_map_Ca.erase(loci_types.first));
		}

		if (type_map_Cb.at(loci_types.second) > 1) {

			type_map_Cb.at(loci_types.second)--;

		} else {

			assert(type_map_Cb.erase(loci_types.second));
		}
	}
}


void TypeContainer::addTypeMap(AncestralTypeMaps & new_type_map) {

	fetched_overlap = false;

	type_map_overlap.clear();
	type_map_overlap.reserve(new_type_map.A.size() + new_type_map.B.size() + new_type_map.C.size());

	for (auto it = new_type_map.A.begin(); it != new_type_map.A.end(); it++) {

		assert(it->second > 0);
		assert (it->first.size() == locus_length_A);

		num_lineages_A += it->second;

		if (!(type_map_A.insert({it->first, it->second}).second)) {

			type_map_A.at(it->first) += it->second;
		} 

		type_map_overlap.push_back(pair <uint,uint> (type_map_A[it->first], it->second)); //old + new, new
	}

	for (auto it = new_type_map.B.begin(); it != new_type_map.B.end(); it++) {

		assert(it->second > 0);
		assert (it->first.size() == locus_length_B);

		num_lineages_B += it->second;

		if (!(type_map_B.insert({it->first, it->second}).second)) {

			type_map_B.at(it->first) += it->second;
		} 

		type_map_overlap.push_back(pair <uint,uint> (type_map_B[it->first], it->second));//old + new, new
	}

	for (auto it = new_type_map.C.begin(); it != new_type_map.C.end(); it++) {

		assert(it->second > 0);
		assert (it->first.size() == haplotype_length);

		num_lineages_C += it->second;

		if (!(type_map_C.insert({it->first, it->second}).second)) {

			type_map_C.at(it->first) += it->second;
		} 

		type_map_overlap.push_back(pair <uint,uint> (type_map_C[it->first], it->second));//old + new, new

 		pair<vector<bool>, vector<bool> > loci_types;
		getLociTypes(it->first, &loci_types);

		if (!(type_map_Ca.insert({loci_types.first, it->second}).second)) {

			type_map_Ca.at(loci_types.first) += it->second;
		} 

		if (!(type_map_Cb.insert({loci_types.second, it->second}).second)) {

			type_map_Cb.at(loci_types.second) += it->second;
		} 
	}

	assert(new_type_map.A.size() + new_type_map.B.size() + new_type_map.C.size() == type_map_overlap.size());
}


vector<pair<uint, uint> > TypeContainer::getTypeMapOverlap() {

	assert(!isEmpty());
	assert(!fetched_overlap);
	fetched_overlap = true;
	
	return type_map_overlap;
}


vector<Haplotype> TypeContainer::getHaplotypesA() {

    vector<Haplotype> return_vector;
    return_vector.reserve(type_map_A.size());

    Haplotype temp_haplotype; 
    temp_haplotype.gamma = "A";

    for (auto it = type_map_A.begin(); it != type_map_A.end(); it++){

        temp_haplotype.type = it->first;
        return_vector.push_back(temp_haplotype);
    }

    return return_vector;
}


vector<Haplotype> TypeContainer::getHaplotypesB() {

	vector<Haplotype> return_vector;
    return_vector.reserve(type_map_B.size());
    uint return_vector_idx = 0; 

    Haplotype temp_haplotype; 
    temp_haplotype.gamma = "B";

    for (auto it = type_map_B.begin(); it != type_map_B.end(); it++){

        temp_haplotype.type = it->first;
        return_vector.push_back(temp_haplotype);
    }

    return return_vector;
}


pair<Haplotype, vector<double> > TypeContainer::sampleType(double theta_A, double theta_B, double rho) {

    vector<SampleInfo> sample_vector;
    sample_vector.reserve(getUniqueCount());

    double norm_const = 0;

    for (auto it = type_map_C.begin(); it != type_map_C.end(); it++) {

        uint count_C = it->second;
        uint count_A = 0;
        uint count_B = 0;

 		pair<vector<bool>, vector<bool> > loci_types;
		getLociTypes(it->first, &loci_types);

        if (type_map_A.count(loci_types.first) > 0) {

        	count_A = type_map_A[loci_types.first];
        }

        if (type_map_B.count(loci_types.second) > 0) {

        	count_B = type_map_B[loci_types.second];
        }

        double probability = (count_C - 1 + (theta_A + theta_B) + rho + count_A + count_B) * count_C;
       	assert (probability > 0);
        
        norm_const += probability;

        SampleInfo sample_info;
        sample_info.probability = probability;
        sample_info.type_p = &it->first;
        sample_info.gamma = "C";
        sample_vector.push_back(sample_info);
    }

    for (auto it = type_map_A.begin(); it != type_map_A.end(); it++) {

        uint count_A = it->second;
		uint count_Ca = 0;

        if (type_map_Ca.count(it->first) > 0) {

        	count_Ca = type_map_Ca[it->first];
        }

        double probability = (count_A - 1 + num_lineages_B + count_Ca + theta_A) * count_A;
        assert (probability > 0);
        
        norm_const += probability;

        SampleInfo sample_info;
        sample_info.probability = probability;
        sample_info.type_p = &it->first;
        sample_info.gamma = "A";
        sample_vector.push_back(sample_info);
    }

    for (TypeMap::iterator it = type_map_B.begin(); it != type_map_B.end(); it++) {

        uint count_B = it->second;
        uint count_Cb = 0;

        if (type_map_Cb.count(it->first) > 0) {

        	count_Cb = type_map_Cb[it->first];
        }

        double probability = (count_B - 1 + num_lineages_A + count_Cb + theta_B) * count_B;
        assert (probability > 0);
        
        norm_const += probability;

        SampleInfo sample_info;
        sample_info.probability = probability;
        sample_info.type_p = &it->first;
        sample_info.gamma = "B";
        sample_vector.push_back(sample_info);
    }

    double threshold = (*uniform_01_sampler_p)(*mt_rng_pt) * norm_const;
    double cumsum = 0;

    for (auto it = sample_vector.begin(); it != sample_vector.end(); it++) {

    	cumsum += it->probability;

    	if (cumsum > threshold) {

    		Haplotype sampled_haplotype;
    		sampled_haplotype.gamma = it->gamma;
    		sampled_haplotype.type = *(it->type_p);

    		vector <double> prob_vec(2,0);
    		prob_vec[0] = it->probability / norm_const;
    		prob_vec[1] = norm_const;

    		return pair<Haplotype, vector<double> > (sampled_haplotype, prob_vec);
    	}
    }

	assert(false);
}


vector<vector<int> > TypeContainer::getDifferences(Haplotype & haplotype) {

	vector <int> haplotype_counts;
	vector <int> difference_A;
	vector <int> difference_B;

	if (haplotype.gamma == "A") {

		haplotype_counts.reserve(type_map_A.size() + type_map_Ca.size() + num_lineages_B);
		difference_A.reserve(type_map_A.size() + type_map_Ca.size() + num_lineages_B);
		difference_B.reserve(type_map_A.size() + type_map_Ca.size() + num_lineages_B);	

		map < uint, uint > redundancy_map;
		calculateDifference(&redundancy_map, haplotype, type_map_A);
		calculateDifference(&redundancy_map, haplotype, type_map_Ca);

		for (auto rit = redundancy_map.rbegin(); rit != redundancy_map.rend(); rit++) {

			haplotype_counts.push_back(rit->second);
			difference_A.push_back(rit->first);
			difference_B.push_back(-1);
		}

		if (num_lineages_B > 0) {

			haplotype_counts.push_back(num_lineages_B);
			difference_A.push_back(-1);
			difference_B.push_back(-1);
		}

	} else if (haplotype.gamma == "B") {

		haplotype_counts.reserve(type_map_B.size() + type_map_Cb.size() + num_lineages_A);
		difference_A.reserve(type_map_B.size() + type_map_Cb.size() + num_lineages_A);
		difference_B.reserve(type_map_B.size() + type_map_Cb.size() + num_lineages_A);	

		map < uint, uint > redundancy_map;
		calculateDifference(&redundancy_map, haplotype, type_map_B);
		calculateDifference(&redundancy_map, haplotype, type_map_Cb);

		for (auto rit = redundancy_map.rbegin(); rit != redundancy_map.rend(); rit++) {

			haplotype_counts.push_back(rit->second);
			difference_A.push_back(-1);
			difference_B.push_back(rit->first);
		}

		if (num_lineages_A > 0) {

			haplotype_counts.push_back(num_lineages_A);
			difference_A.push_back(-1);
			difference_B.push_back(-1);
		}

	} else {

		assert (haplotype.gamma == "C");

		haplotype_counts.reserve(type_map_A.size() + type_map_B.size() + type_map_C.size());
		difference_A.reserve(type_map_A.size() + type_map_B.size() + type_map_C.size());
		difference_B.reserve(type_map_A.size() + type_map_B.size() + type_map_C.size());

		map<pair<int, int>, uint, ComparePair<pair<int, int> > > redundancy_map;

		for (TypeMap::iterator it = type_map_A.begin(); it != type_map_A.end(); it++) {

			uint differences = 0;

			for (uint i = 0; i < it->first.size(); i++) {

				if (haplotype.type[i] != it->first[i]) {

					differences++;
				}
			}

			pair<int, int> diff_key(differences, -1);

			if (!(redundancy_map.insert({diff_key, it->second}).second)) {

				redundancy_map.at(diff_key) += it->second;
			} 
		}

		for (TypeMap::iterator it = type_map_B.begin(); it != type_map_B.end(); it++) {

			uint differences = 0;

			for (uint i = 0; i < it->first.size(); i++) {

				if (haplotype.type[i + second_locus_start_idx] != it->first[i]) {

					differences++;
				}
			}

			pair<int, int> diff_key(-1, differences);

			if (!(redundancy_map.insert({diff_key, it->second}).second)) {

				redundancy_map.at(diff_key) += it->second;
			} 
		}

		for (TypeMap::iterator it = type_map_C.begin(); it != type_map_C.end(); it++) {

			uint differencesA = 0;
			uint differencesB = 0;

			for (uint i = 0; i < it->first.size(); i++) {

				if (haplotype.type[i] != it->first[i]) {

					if (i < second_locus_start_idx) {

						differencesA++;

					} else {

						differencesB++;
					}
				}
			}

			pair<int, int> diff_key(differencesA, differencesB);

			if (!(redundancy_map.insert({diff_key, it->second}).second)) {

				redundancy_map.at(diff_key) += it->second;
			} 
		}

		for (auto rit = redundancy_map.rbegin(); rit != redundancy_map.rend(); rit++) {

			haplotype_counts.push_back(rit->second);
			difference_A.push_back(rit->first.first);
			difference_B.push_back(rit->first.second);
		}
	}

	vector<vector<int> > difference_vector(3);
	difference_vector[0] = haplotype_counts;
	difference_vector[1] = difference_A;
	difference_vector[2] = difference_B;

	return difference_vector;
}


void TypeContainer::calculateDifference(map < uint, uint > * redundancy_map, Haplotype & haplotype, TypeMap & typemap) {

	for (auto it = typemap.begin(); it != typemap.end(); it++) {

		uint differences = 0;

		for (uint i = 0; i < it->first.size(); i++) {

			if (haplotype.type[i] != it->first[i]) {

				differences++;
			}
		}

		if (!(redundancy_map->insert({differences, it->second}).second)) {

			redundancy_map->at(differences) += it->second;
		}
	}	
}





