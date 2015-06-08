
/*
Sample.hpp - This file is part of the Coalescenator (v1.0.0)


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


#ifndef __coalescenator__Samples_hpp
#define __coalescenator__Samples_hpp

#include <string>
#include <vector>
#include <unordered_map>

#include "Utils.hpp"

using namespace std;

typedef vector < vector < vector < uint > > >  SampleContainer;
typedef unordered_map<vector<uint>, uint, VectorHash<vector<uint> > > FullTypeMap;

class Sample {

  protected:

    // Data structures
    vector < AncestralTypeMaps > type_maps;
    vector < double > times;
    vector < double > viral_loads;
    pair<uint, uint> seq_lens;

  public: 

    Sample();
    Sample(vector < AncestralTypeMaps >, vector < double >, vector < double >, pair<uint, uint>);

    // Get functions
    vector < AncestralTypeMaps > getTypeMaps();
    vector < double > getTimes();
    vector < double > getViralLoads();
    pair<uint, uint> getSeqLens();

};

class CompleteSample : public Sample {

  private:

    SampleContainer samples; 
    
    vector <bool> has_missing;
    vector <bool> has_indel;
    vector <bool> has_ambiguous;
    vector <uint> poly_count;

    vector <pair <uint, uint> > loci_coordinates;
    vector <pair <uint, uint> > original_loci_coordinates;

    vector < FullTypeMap > full_type_maps;
    // vector < TypeMap > bi_allelic_type_maps;

    // Generic function for filtering
    template <typename T>
    uint genericFilterFunction(vector<T>, uint, uint);

  public:

    // Constructor
    CompleteSample(string, string, string);
    
    void parseBreakpoints(uint);
    void parseBreakpoints(string);

    void filterAmbiguousSites();

    // Filter out all missing data from sample 
    void filterMissingData();

    // Filter out all indels 
    void filterIndels();

    void filterNonBiAllelic();

    pair<Sample, vector <uint> > getSubSample(pair<uint,uint> , pair<uint,uint>);

    vector <pair <uint, uint> > getOriginalLociCoordinates();
    vector <pair <uint, uint> > getLociCoordinates();
};


#endif
