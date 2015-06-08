
/*
main.cpp - This file is part of the Coalescenator (v1.0.0)


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
#include <fstream>
#include <sstream>
#include <string>

#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>

#include "Utils.hpp"
#include "Sample.hpp"
#include "SamplerOutput.hpp"
#include "ImportanceSampler.hpp"
#include "Matrix3d.hpp"

namespace po = boost::program_options;

using namespace std;

int main (int argc, char * const argv[]) {
       
    OptionsContainer options_variables;

    try {

        po::options_description general("== General options ==", 160);
        general.add_options()
            ("help,h", "produce help message")      
            ("num-threads,p", po::value<uint>(&options_variables.num_threads)->default_value(1), "number of importance sampling threads")
            ("seed,s", po::value<uint>(&options_variables.seed)->default_value(time(NULL), "unix time"), "seed for pseudo-random number generator") 
            ("output-file-prefix", po::value<string>(&options_variables.output_file_prefix)->default_value("coalescenator"), "output file prefix") 
		;

        po::options_description input("== Input options ==", 160);
        input.add_options()
            ("fasta-list,f", po::value<string>(&options_variables.fasta_list)->required(), "list of serial samples in fasta format from last to earliest sample (sample1.fa,sample2.fa,...,sampleN.fa)")          
            ("time-list,t", po::value<string>(&options_variables.time_list)->required(), "list of sample times in the same order as samples (time1,time2,...,timeN)") 
            ("viral-list,v", po::value<std::string>(&options_variables.viral_list)->required(), "list of viral loads in the same order as samples (viral1,viral2,...,viralN)") 
            ("generation-time,g", po::value<double>(&options_variables.generation_time)->required(), "generation time")          
            ("coordinate-list", po::value<string>(&options_variables.coordinate_list), "list of 1-based loci coordinates ([l1_left-l1_right],...,[ln_left-ln_right]) or list of breakpoints - 1-based position after breakpoint (bp1,...bpn). Mutually exclusive with loci-size") 
            ("loci-size", po::value<uint>(&options_variables.loci_size)->default_value(1), "loci size (all sites). Mutually exclusive with coordinate-list") 

        ;

		po::options_description importance("== Importance sampler options ==", 160);
        importance.add_options()            
            ("mc-samples,i", po::value<uint>(&options_variables.mc_samples)->default_value(500000), "number of monto carlo samples")
            ("mutation-list", po::value<string>(&options_variables.mutation_list)->default_value("1e-6,5e-6,1e-5,2.5e-5,5e-5,7.5e-5,1e-4"), "list of comma seperated mutation rates (mutations per generation per site)")   
            ("mutation-driving", po::value<double>(&options_variables.mutation_driving)->default_value(2.5e-5,"2.5e-5"), "driving value for the mutation rate (mutations per generation per site)")
            ("scale-list", po::value<string>(&options_variables.scale_list)->default_value("250,500,750,1000,2000,4000"), "list of comma seperated scaling parameters (proportionality constant bewteen viral load and effective population size)")   
            ("scale-driving", po::value<double>(&options_variables.scale_driving)->default_value(1000), "driving value for the scaling parameter (proportionality constant bewteen viral load and effective population size)")
            ("rho-list", po::value<string>(&options_variables.rho_list)->default_value("1e-10,1e-9,1e-8,1e-7,1e-6,5e-6,1e-5,5e-5,1e-4,5e-4,1e-3"), "list of comma seperated recombination rates (recombinations per generation per site)")   
            ("rho-driving", po::value<double>(&options_variables.rho_driving)->default_value(5e-5,"5e-5"), "driving value for the recombination rate (recombinations per generation per site)")
        ;    
        
        po::options_description composite("== Composite likelihood options ==", 160);
        composite.add_options()            
            ("min-comp-dist", po::value<uint>(&options_variables.min_dist)->default_value(1), "minimum distance (number of nucleotides) between loci in the composite likelihood estimation")
            ("max-comp-dist", po::value<uint>(&options_variables.max_dist)->default_value(400), "maximum distance (number of nucleotides) between loci in the composite likelihood estimation")
        ; 

        po::options_description desc("\n#### Coalescenator ####");
        desc.add(general).add(input).add(importance).add(composite);
                
        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        
        if (vm.count("help") || argc == 1) {
            cout << desc << endl;
            return 1;
        }
        
        po::notify(vm);

        options_variables.qadrature_points = 16;

        cout << "\n[" << getLocalTime() << "] Seed used for pseudo-random generator: " << options_variables.seed << "\n" << endl;

        CompleteSample completeSample(options_variables.fasta_list, options_variables.time_list, options_variables.viral_list);
        
        if (vm.count("coordinate-list")) {

            completeSample.parseBreakpoints(options_variables.coordinate_list);      
        
        } else {

            assert (!vm.count("coordinate-list"));
            completeSample.parseBreakpoints(options_variables.loci_size);      
        }    

        completeSample.filterAmbiguousSites();
        completeSample.filterIndels();
        completeSample.filterNonBiAllelic();

        vector < pair <uint,uint > > loci_coordinates = completeSample.getLociCoordinates();


        if (loci_coordinates.size() == 0) {

           cout << "\n[" << getLocalTime() << "] All loci was removed by the filtering process - exiting!" << endl; 
           exit(0);
        }

        cout << "\n[" << getLocalTime() << "] Number of loci: " << loci_coordinates.size() << endl;
        cout << "[" << getLocalTime() << "] List of loci lengths: \n" << endl;

        cout << "\t" << loci_coordinates[0].second - loci_coordinates[0].first + 1;

        for (uint i = 1; i < loci_coordinates.size(); i++) {

            cout << ", " << loci_coordinates[i].second - loci_coordinates[i].first + 1;

        }

        cout << endl;

        assert (options_variables.mutation_driving > 0);
        assert (options_variables.scale_driving > 0);
        assert (options_variables.rho_driving > 0);
        assert (options_variables.generation_time > 0);
        
        ModelParameters model_parameters;

        model_parameters.driving_mu = options_variables.mutation_driving;
        model_parameters.driving_lambda = options_variables.scale_driving;
        model_parameters.driving_rho = options_variables.rho_driving;
        model_parameters.tau = options_variables.generation_time;

        vector <string> mu_list_string;  
        vector <string> lambda_list_string;
        vector <string> rho_list_string;

        boost::split(mu_list_string, options_variables.mutation_list, boost::is_any_of(","));
        boost::split(lambda_list_string, options_variables.scale_list, boost::is_any_of(","));
        boost::split(rho_list_string, options_variables.rho_list, boost::is_any_of(","));

        vector <double> mu_list(mu_list_string.size(), 0);  
        vector <double> lambda_list(lambda_list_string.size() ,0);
        vector <double> rho_list(rho_list_string.size() ,0);

        for (int i = 0; i < mu_list_string.size(); i++) {

            mu_list[i] = ::atof(mu_list_string[i].c_str());
            assert (mu_list[i] >= 0);
        }

        for (int i = 0; i < lambda_list_string.size(); i++) {

            lambda_list[i] = ::atof(lambda_list_string[i].c_str());
            assert (lambda_list[i] > 0);
        }

        for (int i = 0; i < rho_list_string.size(); i++) {

            rho_list[i] = ::atof(rho_list_string[i].c_str());
            assert (rho_list[i] >= 0);
        } 

        cout << "\n[" << getLocalTime() << "] Parameter list: " << endl;
        cout << "\n\t Mutation rate(s): " << mu_list << endl;
        cout << "\t Driving mutation rate: " << options_variables.mutation_driving << endl;
        cout << "\n\t Recombination rate(s): " << rho_list << endl;
        cout << "\t Driving recombination rate: " << options_variables.rho_driving << endl;
        cout << "\n\t Scaling parameter(s): " << lambda_list << endl;
        cout << "\t Driving scaling parameter: " << options_variables.scale_driving << endl;

        model_parameters.target_mu = mu_list;
        model_parameters.target_lambda = lambda_list;
        model_parameters.target_rho = rho_list;

        vector <pair <uint, uint> > original_loci_coordinates = completeSample.getOriginalLociCoordinates();
        assert (original_loci_coordinates.size() == loci_coordinates.size());

        if (loci_coordinates.size() == 1) {

           cout << "\n[" << getLocalTime() << "] Only a single loci left no remcobination estimation possible - exiting!" << endl; 
           exit(0);

        }

        uint locus_counter = 0;

        for (auto it_first = original_loci_coordinates.begin(); it_first != original_loci_coordinates.end(); it_first++) {

            for (auto it_second = it_first + 1; it_second != original_loci_coordinates.end(); it_second++) {

                uint loci_distance = it_second->first - it_first->second - 1;
                assert (loci_distance >= 0);

                if ((loci_distance < options_variables.min_dist) or (loci_distance > options_variables.max_dist)) {

                    continue;
                }

                locus_counter++;
            }
        }

        cout << "\n[" << getLocalTime() << "] For each of the " << locus_counter << " loci combination(s) generate " << options_variables.mc_samples << " importance samples. " << options_variables.mc_samples*locus_counter << " samples in total.\n" << endl;

        Matrix3d<double> composite_likelihood(model_parameters.target_mu.size(), model_parameters.target_rho.size(), model_parameters.target_lambda.size(), 0);
        Matrix3d<double> mean_tmrca(model_parameters.target_mu.size(), model_parameters.target_rho.size(), model_parameters.target_lambda.size(), 0);

        stringstream output_file_pairwise_ss; 
        output_file_pairwise_ss << options_variables.output_file_prefix << "_pairwise_loglikelihood.txt";
        string output_file_pairwise_str = output_file_pairwise_ss.str();

        ofstream output_file_pairwise;
        output_file_pairwise.open(output_file_pairwise_str.c_str());
        assert (output_file_pairwise.is_open());
        
        output_file_pairwise << "locus1\tlocus2\tmu\trho\tlambda\tloglikelihood\ttmrca\tlineages_remaining\tloglikelihoodSq\n";
        
        for (auto i = 0; i != loci_coordinates.size(); i++) {

            for (auto j = i + 1; j != loci_coordinates.size(); j++) {

                uint loci_distance = original_loci_coordinates[j].first - original_loci_coordinates[i].second - 1;
                assert (loci_distance >= 0);

                if ((loci_distance < options_variables.min_dist) or (loci_distance > options_variables.max_dist)) {

                    continue;
                }

                pair <Sample, vector<uint> > sampleContainer = completeSample.getSubSample(loci_coordinates[i], loci_coordinates[j]);

                Sample subSample = sampleContainer.first;
                uint num_times = subSample.getTimes().size();
                
                //change target and driving rho to be proportional to loci
                model_parameters.target_rho = rho_list;
                for (auto it_rho = model_parameters.target_rho.begin(); it_rho != model_parameters.target_rho.end(); it_rho++) {

                    *it_rho *= loci_distance + 1;
                    
                    // if (*it_rho > 0.5) {

                    //     *it_rho = 0.5;
                    // }
                } 
                model_parameters.driving_rho = options_variables.rho_driving;
                model_parameters.driving_rho *= loci_distance + 1;

                locus_counter--;
                stringstream locus1_string;
                locus1_string << "(" << original_loci_coordinates[i].first << ", " << original_loci_coordinates[i].second << ")";
                stringstream locus2_string;
                locus2_string << "(" << original_loci_coordinates[j].first << ", " << original_loci_coordinates[j].second << ")";

                cout << "[" << getLocalTime() << "] Working on locus " << locus1_string.str() << " and " << locus2_string.str() << " having " << sampleContainer.second[0] << ", " << sampleContainer.second[1] << " and " << sampleContainer.second[2] << " sequences with ancentry A, B and C respectively taken at " << num_times << " different time points. " << locus_counter << " remaining pairwise combinations." << endl;

                ImportanceSampler importance_sampler;

                SamplerOutput<double> results = importance_sampler.runSampler(&subSample, model_parameters, options_variables.mc_samples, options_variables.num_threads, options_variables.qadrature_points, options_variables.seed);

                composite_likelihood.addMatrixElementWise(results.likelihood());
                mean_tmrca.logAddMatrixElementWise(results.tmrca());

                for (uint i = 0; i < model_parameters.target_mu.size(); i++) {

                    for (uint j = 0; j < model_parameters.target_rho.size(); j++) {

                        for (uint k = 0; k < model_parameters.target_lambda.size(); k++) {

                            output_file_pairwise << locus1_string.str() << "\t" << locus2_string.str() << "\t" << model_parameters.target_mu[i] << "\t" << rho_list[j] << "\t" << model_parameters.target_lambda[k] << "\t" << results.likelihood().getElement(i,j,k) << "\t" << exp(results.tmrca().getElement(i,j,k)) << "\t";
                            
                            if (results.lineages_remaining().size() > 0) {
                        
                                output_file_pairwise << exp(results.lineages_remaining()[0].getElement(i,j,k));

                                for (auto it_lin = results.lineages_remaining().begin() + 1; it_lin != results.lineages_remaining().end(); it_lin++) {

                                    output_file_pairwise << "," << exp(it_lin->getElement(i,j,k));
                                }
                            }

                            output_file_pairwise << "\t" << results.likelihoodSq().getElement(i,j,k);
                            output_file_pairwise << "\n";
                        }
                    }
                }
            }
        }

        output_file_pairwise.close();
        assert (locus_counter == 0);

        stringstream output_file_composite_ss; 
        output_file_composite_ss << options_variables.output_file_prefix << "_composite_loglikelihood.txt";
        string output_file_composite_str = output_file_composite_ss.str();

        ofstream output_file_composite;
        output_file_composite.open(output_file_composite_str.c_str());
        assert (output_file_composite.is_open());

        mean_tmrca.subtractConstant(log(loci_coordinates.size() * (loci_coordinates.size() - 1)/2));

        output_file_composite << "mu\trho\tlambda\tlikelihood\ttmrca\n";

        for (uint i = 0; i < model_parameters.target_mu.size(); i++) {

            for (uint j = 0; j < model_parameters.target_rho.size(); j++) {

                for (uint k = 0; k < model_parameters.target_lambda.size(); k++) {
                             
                    output_file_composite << model_parameters.target_mu[i] << "\t" << rho_list[j] << "\t" << model_parameters.target_lambda[k] << "\t" << composite_likelihood.getElement(i,j,k) << "\t" << exp(mean_tmrca.getElement(i,j,k)) << "\n";              
                }
            }
        }

        output_file_composite.close();
	}
    
    catch(std::exception& e) {
        cerr << "ERROR: " << e.what() << endl;
        return 1;
    }
    
    catch(...) {
        cerr << "Exception of unknown type!" << endl;
        return 1;
    }

    cout << "[" << getLocalTime() << "] Coalescenated succesfully." << endl;;
    return 0;   
}