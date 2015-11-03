
/*
Sample.cpp - This file is part of the Coalescenator (v1.0.0)


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


#include "Sample.hpp"

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>

#include <boost/algorithm/string.hpp>

#include "Utils.hpp"

using namespace std;

Sample::Sample(){}

Sample::Sample(vector < AncestralTypeMaps > type_maps_in, vector < double > times_in, vector < double > viral_loads_in, pair<uint, uint> seq_lens_in) {

	type_maps = type_maps_in;
	times = times_in;
	viral_loads = viral_loads_in;
	seq_lens = seq_lens_in;
}


vector < AncestralTypeMaps > Sample::getTypeMaps() {

	return type_maps;
}


vector < double > Sample::getTimes() {

	return times;
}


vector < double > Sample::getViralLoads() {

	return viral_loads;
}


pair<uint, uint> Sample::getSeqLens() {

	return seq_lens;
}


// Constructor
CompleteSample::CompleteSample (string fasta_in, string time_in, string viral_in) {

	// Read fasta filenames
	vector <string> sample_string_vec;
	vector <string> time_string_vec;
    vector <string> viral_string_vec;

	boost::split(sample_string_vec, fasta_in, boost::is_any_of(","));
	boost::split(time_string_vec, time_in, boost::is_any_of(","));
	boost::split(viral_string_vec, viral_in, boost::is_any_of(","));

	//vector <string> viral_string_vec(time_string_vec.size(),"1");

	assert (time_string_vec.size() == sample_string_vec.size());
	assert (viral_string_vec.size() == sample_string_vec.size());

	vector <string> sorted_sample_vec;
	sorted_sample_vec.reserve(time_string_vec.size());

	vector <double> sorted_time_vec;
	sorted_time_vec.reserve(time_string_vec.size());

	vector <double> sorted_viral_vec;
	sorted_viral_vec.reserve(time_string_vec.size());

	// Sort the sample, time and viral lists according to sample times (bubble sort)
	for (auto i = 0; i < time_string_vec.size(); i++) {

		if (sorted_time_vec.size() == 0) {

			sorted_sample_vec.push_back(sample_string_vec[i]);
			sorted_time_vec.push_back(::atof(time_string_vec[i].c_str()));
			sorted_viral_vec.push_back(::atof(viral_string_vec[i].c_str()));

		} else {

			sorted_sample_vec.push_back(sample_string_vec[i]);
			sorted_time_vec.push_back(::atof(time_string_vec[i].c_str()));
			sorted_viral_vec.push_back(::atof(viral_string_vec[i].c_str()));

			for (auto j = sorted_time_vec.size() - 1; j > 0; j--) {

				if (sorted_time_vec[j-1] > sorted_time_vec[j]) {

					string temp_sample = sorted_sample_vec[j-1];
					sorted_sample_vec[j-1] = sorted_sample_vec[j];
					sorted_sample_vec[j] = temp_sample;

					double temp_time = sorted_time_vec[j-1];
					sorted_time_vec[j-1] = sorted_time_vec[j];
					sorted_time_vec[j] = temp_time;

					double temp_viral = sorted_viral_vec[j-1];
					sorted_viral_vec[j-1] = sorted_viral_vec[j];
					sorted_viral_vec[j] = temp_viral;
				}
			}
		}
	}

	times = sorted_time_vec;
	viral_loads = sorted_viral_vec;

	assert (times.size() == viral_loads.size());
	assert (times.size() == sorted_sample_vec.size());

	SampleContainer temp_sample_list;
	vector < unordered_map <string, uint> > type_count_map;
	unordered_set<vector<uint>, VectorHash<uint> > unique_set;

	string fasta_id_char = ">";

	uint read_counter = 0;
  	uint length_test = 0;

	// Loop over fasta files
	for (auto i = 0; i < sorted_sample_vec.size(); i++) {

		string line;
  		ifstream file_handle(sorted_sample_vec[i].c_str());
  		assert (file_handle.is_open());

  		FullTypeMap sample_type_maps;

  		cout << "[" << getLocalTime() << "] Parsing: " << sorted_sample_vec[i] << " taken at time " << times[i] << " and with a measured concentration of " << viral_loads[i] << "." << endl;

		while (file_handle.good()) {

  			getline (file_handle, line);

  			if ((line.size() > 0) and (fasta_id_char[0] != line[0])) {

  				if (read_counter == 0) {

	  				 length_test = line.size();

	  			} else {

	  				assert (length_test == line.size());
	  			}

  				read_counter += 1;

  				bool left_missing = true;

  				if (type_count_map.size() == 0) {

      				for (auto j = 0; j < line.size(); j++) {

      					unordered_map <string, uint> temp_map;
      					type_count_map.push_back(temp_map);

      					has_missing.push_back(false);
      					has_indel.push_back(false);
      					has_ambiguous.push_back(false);
      				}
      			}

  				vector < uint > cur_sequence;
  				cur_sequence.reserve(has_indel.size());
      			vector <bool> has_indel_cur_seq(has_indel.size(), false);

  				for (auto j = 0; j < line.size(); j++) {

  					if (line.substr(j,1) == "-" and left_missing == true) {

  						has_missing[j] = true;
  						cur_sequence.push_back(5);
  						continue;
  					}

  					left_missing = false;

  					if (line.substr(j,1) == "-") {

  						has_indel_cur_seq[j] = true;
  						cur_sequence.push_back(4);

  					} else if (line.substr(j,1) == "N") {

  						has_ambiguous[j] = true;
  						cur_sequence.push_back(6);

  					} else {

  						if (type_count_map[j].count(line.substr(j,1)) > 0) {

      						cur_sequence.push_back(type_count_map[j][line.substr(j,1)]);

      					} else {

      						int num_types = type_count_map[j].size();
      						assert(type_count_map[j].insert({line.substr(j,1), num_types}).second);

      						cur_sequence.push_back(type_count_map[j][line.substr(j,1)]);
      					}
  					}
      			}

  				for (int j = line.size() - 1; j > -1; j--) {

  					if (line.substr(j,1) == "-") {

  						has_missing[j] = true;
  						has_indel_cur_seq[j] = false;
  						cur_sequence[j] = 5;
  						continue;
  					}

  					break;
  				}

      			assert (has_missing.size() == cur_sequence.size());
      			assert (has_indel.size() == cur_sequence.size());
      			assert (type_count_map.size() == cur_sequence.size());

      			for (auto j = 0; j < has_indel_cur_seq.size(); j++) {

      				if (has_indel_cur_seq[j]) {

      					has_indel[j] = true;
      				}
      			}

				unique_set.insert(cur_sequence);

      			if (!(sample_type_maps.insert({cur_sequence, 1}).second)) {

      				sample_type_maps.at(cur_sequence) += 1;
      			}
  			}
    	}

    	full_type_maps.push_back(sample_type_maps);
    	file_handle.close();
	}

	for (int i = 0; i < has_missing.size(); i++) {

		if (has_indel[i] == true) {

			poly_count.push_back(type_count_map[i].size() + 1);

		} else {

			poly_count.push_back(type_count_map[i].size());
		}
	}

	assert (has_ambiguous.size() == has_indel.size());
	assert (poly_count.size() == has_indel.size());
	assert (has_missing.size() == has_indel.size());

	cout << "\n[" << getLocalTime() << "] " << read_counter << " read parsed in total for " << times.size() << " serially sampled samples." << endl;
	cout << "[" << getLocalTime() << "] " << has_missing.size() << " sites in the supplied chromosomes with " << unique_set.size() << " unique types in total." << endl;
}


void CompleteSample::parseBreakpoints(string input_string) {

    assert (loci_coordinates.size() == 0);
    assert (original_loci_coordinates.size() == 0);

    vector <pair <uint, uint> > temp_loci_coord;

    vector <string> input_check_left_bracket;
	boost::split(input_check_left_bracket, input_string, boost::is_any_of("["));

	if (input_check_left_bracket.size() > 1) {

    	vector <string> input_check_right_bracket;
		boost::split(input_check_right_bracket, input_string, boost::is_any_of("]"));
		assert (input_check_left_bracket.size() == input_check_right_bracket.size());

	    vector <string> loci_string_split;
	    boost::split(loci_string_split, input_string, boost::is_any_of(","));
	    assert (loci_string_split.size() > 0);
	    assert (loci_string_split.size()== input_check_left_bracket.size() - 1);

	    for (auto it = loci_string_split.begin(); it != loci_string_split.end(); it++) {

	        assert (it->size() > 0);
	        assert (it->substr(0,1)  == "[");
	        assert (it->substr(it->size() - 1,1)  == "]");

	        it->erase(it->begin());
	        it->erase(it->end() - 1);

	        vector <string> loci_coord_split;
	        boost::split(loci_coord_split, *it, boost::is_any_of("-"));
	        assert (loci_coord_split.size() == 2);

	        pair<uint, uint> temp_pair(::atof(loci_coord_split[0].c_str()), ::atof(loci_coord_split[1].c_str()));
	        assert (temp_pair.first <= temp_pair.second);

	        if (temp_loci_coord.size() > 0) {

	            assert (temp_loci_coord.back().second < temp_pair.first);

	        }

	        temp_loci_coord.push_back(temp_pair);
	    }

	} else {

	    vector <string> break_string_split;
	    boost::split(break_string_split, input_string, boost::is_any_of(","));
	    assert (break_string_split.size() > 0);

	    vector <uint> break_points;

	    for (uint i = 0; i < break_string_split.size(); i++) {

	        break_points.push_back(::atof(break_string_split[i].c_str()));
	        assert (break_points.back() > 1);
	        assert (break_points.back() <= has_missing.size());

	        if (i == 0) {

	            pair<uint, uint> temp_pair(1, break_points[i] - 1);
	            assert (temp_pair.first <= temp_pair.second);
	            temp_loci_coord.push_back(temp_pair);

	        } else {

	            pair<uint, uint> temp_pair(break_points[i - 1], break_points[i] - 1);
	            assert (temp_pair.first <= temp_pair.second);
	            temp_loci_coord.push_back(temp_pair);
	        }
	    }

	    pair<uint, uint> temp_pair(break_points.back(), has_missing.size());
	    assert (temp_pair.first <= temp_pair.second);
	    temp_loci_coord.push_back(temp_pair);
	}

    assert (temp_loci_coord.front().first > 0);
	assert (temp_loci_coord.back().second <= has_missing.size());

    loci_coordinates = temp_loci_coord;
    original_loci_coordinates = temp_loci_coord;
}


void CompleteSample::parseBreakpoints(uint loci_size) {

    assert (loci_coordinates.size() == 0);
    assert (original_loci_coordinates.size() == 0);
    assert (loci_size <= has_missing.size());

    vector <pair <uint, uint> > temp_loci_coord;

    uint left_most = 1;
    uint right_most = loci_size;

    while (right_most <= has_missing.size()) {

        temp_loci_coord.push_back(pair <uint, uint> (left_most, right_most));
        left_most = right_most + 1;
        right_most = left_most + loci_size - 1;
    }

    if (temp_loci_coord.back().second < has_missing.size()) {

    	temp_loci_coord.push_back(pair <uint, uint> (temp_loci_coord.back().second + 1, has_missing.size()));
    }

    loci_coordinates = temp_loci_coord;
    original_loci_coordinates = temp_loci_coord;
}


vector <pair <uint, uint> > CompleteSample::getOriginalLociCoordinates() {

	return original_loci_coordinates;
}


vector <pair <uint, uint> > CompleteSample::getLociCoordinates() {

	return loci_coordinates;
}


void CompleteSample::filterAmbiguousSites() {

	uint unique_count = genericFilterFunction(has_ambiguous, 0, 0);
	cout << "[" << getLocalTime() << "] " << has_ambiguous.size() << " sites and " << unique_count << " unique types remaining after filtering ambiguous sites." << endl;
}


void CompleteSample::filterMissingData() {

	assert(false);

	uint unique_count = genericFilterFunction(has_missing, 0, 0);
	cout << "[" << getLocalTime() << "] " << has_missing.size() << " sites and " << unique_count << " unique types remaining after filtering missing data." << endl;
}


void CompleteSample::filterIndels() {

	uint unique_count = genericFilterFunction(has_indel, 0, 0);
	cout << "[" << getLocalTime() << "] " << has_indel.size() << " sites and " << unique_count << " unique types remaining after filtering indels." << endl;
}


void CompleteSample::filterNonBiAllelic() {

	uint unique_count = genericFilterFunction(poly_count, 2, 2);
	cout << "[" << getLocalTime() << "] " << poly_count.size() << " sites and " << unique_count << " unique types remaining after filtering non-bi-allelic positions." << endl;
}

// Generic function for filtering
template <typename T>
uint CompleteSample::genericFilterFunction(vector <T> filter_vec, uint min, uint max) {

	bool first_seq = true;

	vector < bool > has_ambiguous_temp;
	vector < bool > has_missing_temp;
	vector < bool > has_indel_temp;
	vector < uint > poly_count_temp;

	has_ambiguous_temp.reserve(has_ambiguous.size());
	has_missing_temp.reserve(has_missing.size());
	has_indel_temp.reserve(has_indel.size());
	poly_count_temp.reserve(poly_count.size());

	unordered_set<vector<uint>, VectorHash<uint> > unique_set;
	vector <uint> loci_coordinates_permutation(filter_vec.size(), 0);

	FullTypeMap temp_type_map;

	for (auto it_samp = full_type_maps.begin(); it_samp != full_type_maps.end(); it_samp++) {

		temp_type_map.clear();

		for (auto it_seq = it_samp->begin(); it_seq != it_samp->end(); it_seq++) {

			vector <uint> temp_type;
			temp_type.reserve(it_seq->first.size());

			assert (it_seq->first.size() == filter_vec.size());
			uint cumsum_permutation = 0;

			for (int k = 0; k < filter_vec.size(); k++) {

				if (filter_vec[k] >= min and filter_vec[k] <= max) {

					temp_type.push_back(it_seq->first[k]);

					if (first_seq) {

						has_ambiguous_temp.push_back(has_ambiguous[k]);
						has_missing_temp.push_back(has_missing[k]);
						has_indel_temp.push_back(has_indel[k]);
						poly_count_temp.push_back(poly_count[k]);
						cumsum_permutation += 1;
					}
				}

				if (first_seq) {

					loci_coordinates_permutation[k] = cumsum_permutation;
				}
			}

			first_seq = false;
  			unique_set.insert(temp_type);

		   	if (!(temp_type_map.insert({temp_type, it_seq->second}).second)) {

  				temp_type_map.at(temp_type) += it_seq->second;
  			}
		}

		*it_samp = temp_type_map;
	}

    vector <uint> erase_loci;

	for (auto it = loci_coordinates.begin(); it != loci_coordinates.end(); it++) {

		if (filter_vec[it->first - 1] >= min and filter_vec[it->first - 1] <= max) {

			it->first = loci_coordinates_permutation[it->first - 1];

		} else {

			it->first = loci_coordinates_permutation[it->first - 1] + 1;
		}

		it->second = loci_coordinates_permutation[it->second - 1];

		if (it->first > it->second) {

			assert (it->first == it->second + 1);
			erase_loci.push_back(it - loci_coordinates.begin());
		}
	}

    sort(erase_loci.begin(), erase_loci.end());

    vector <pair <uint, uint> > erased_loci_coordinates;
    erased_loci_coordinates.reserve(erase_loci.size());

    assert (erase_loci.size() <= loci_coordinates.size());
    assert (original_loci_coordinates.size() == loci_coordinates.size());

	for (auto rit = erase_loci.rbegin(); rit != erase_loci.rend(); rit++) {

		erased_loci_coordinates.push_back(pair<uint, uint> (original_loci_coordinates[*rit].first, original_loci_coordinates[*rit].second));
		loci_coordinates.erase(loci_coordinates.begin() + *rit);
		original_loci_coordinates.erase(original_loci_coordinates.begin() + *rit);
	}

	if (!erased_loci_coordinates.empty()) {

		cout << "\n";
	}

	for (auto rit = erased_loci_coordinates.rbegin(); rit != erased_loci_coordinates.rend(); rit++) {

		cout << "[" << getLocalTime() << "] WARNING: Removed loci " << "(" << rit->first << ", " << rit->second << ") due to a lack of relevant sites." << endl;
	}

	if (!erased_loci_coordinates.empty()) {

		cout << "\n";
	}

	assert (has_ambiguous_temp.size() == has_missing_temp.size());
	assert (has_ambiguous_temp.size() == has_indel_temp.size());
	assert (has_ambiguous_temp.size() == poly_count_temp.size());

	has_ambiguous = has_ambiguous_temp;
	has_missing = has_missing_temp;
	has_indel = has_indel_temp;
	poly_count = poly_count_temp;

	return unique_set.size();
}


pair<Sample, vector <uint> > CompleteSample::getSubSample(pair<uint,uint> first_locus_idx, pair<uint,uint> second_locus_idx) {

	vector <uint> type_counts(3,0);

	vector < AncestralTypeMaps > temp_type_maps;
	vector<double> temp_times;
	vector<double> temp_viral_loads;

	for (auto i = 0; i < full_type_maps.size(); i++) {

		AncestralTypeMaps temp_type_map;
		uint seq_counter = 0;

		for (auto it_seq = full_type_maps.at(i).begin(); it_seq != full_type_maps.at(i).end(); it_seq++) {

			bool ancenstralA = true;
			bool ancenstralB = true;

			for (uint j = first_locus_idx.first - 1; j < first_locus_idx.second; j++) {

				if (it_seq->first[j] == 5) {

					ancenstralA = false;
					break;
				}
			}

			for (uint j = second_locus_idx.first - 1; j < second_locus_idx.second; j++) {

				if (it_seq->first[j] == 5) {

					ancenstralB = false;
					break;
				}
			}

			if (ancenstralA == false and ancenstralB == false) {

				continue;
			}

			if (ancenstralA == true and ancenstralB == false) {

				vector <bool> temp_type(it_seq->first.begin() + first_locus_idx.first - 1, it_seq->first.begin() + first_locus_idx.second);

				if (!(temp_type_map.A.insert({temp_type, it_seq->second}).second)) {

					temp_type_map.A.at(temp_type) += it_seq->second;
				}

				type_counts[0] += it_seq->second;
				seq_counter++;
			}

			if (ancenstralA == false and ancenstralB == true) {

				vector <bool> temp_type(it_seq->first.begin() + second_locus_idx.first - 1, it_seq->first.begin() + second_locus_idx.second);

				if (!(temp_type_map.B.insert({temp_type, it_seq->second}).second)) {

					temp_type_map.B.at(temp_type) += it_seq->second;
				}

				type_counts[1] += it_seq->second;
				seq_counter++;
			}

			if (ancenstralA == true and ancenstralB == true) {

				vector <bool> temp_type(it_seq->first.begin() + first_locus_idx.first - 1, it_seq->first.begin() + first_locus_idx.second);
				temp_type.insert(temp_type.end(), it_seq->first.begin() + second_locus_idx.first - 1, it_seq->first.begin() + second_locus_idx.second);

				if (!(temp_type_map.C.insert({temp_type, it_seq->second}).second)) {

					temp_type_map.C.at(temp_type) += it_seq->second;
				}

				type_counts[2] += it_seq->second;
				seq_counter++;
			}
		}

		if (seq_counter != 0) {

			temp_type_maps.push_back(temp_type_map);
			temp_times.push_back(times.at(i));
			temp_viral_loads.push_back(viral_loads.at(i));
		}
	}

	if (temp_type_maps.empty()) {

		cout << "ERROR: All sequences in both loci have missing (increase --loci-size or specify another region)" << endl;
		exit(1);
	}


	assert (type_counts[0] > 0 or type_counts[1] > 0 or type_counts[2] > 0);

	pair <uint, uint> seq_lengths(first_locus_idx.second - first_locus_idx.first + 1, second_locus_idx.second - second_locus_idx.first + 1);
	assert (seq_lengths.first > 0 and seq_lengths.second > 0);

	Sample subSample(temp_type_maps, temp_times, temp_viral_loads, seq_lengths);

	return pair<Sample, vector<uint> > (subSample, type_counts);
}
