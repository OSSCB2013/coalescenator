
/*
TypeContainer.hpp - This file is part of the Coalescenator (v1.0.0)


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


#ifndef __coalescenator__TypeContainer_hpp
#define __coalescenator__TypeContainer_hpp

#include <vector>
#include <random>

#include <boost/functional/hash.hpp>

#include "Utils.hpp"

using namespace std;

typedef uniform_real_distribution<double> uniform_01_sampler_t;

class TypeContainer {

  private:

    mt19937* mt_rng_pt;
    uniform_01_sampler_t* uniform_01_sampler_p;

    TypeMap type_map_A;
    TypeMap type_map_B;
    TypeMap type_map_C;
    TypeMap type_map_Ca;
    TypeMap type_map_Cb;

    uint num_lineages_A;
    uint num_lineages_B;
    uint num_lineages_C;

    uint second_locus_start_idx;
    uint haplotype_length;
    uint locus_length_A;
    uint locus_length_B;

    bool fetched_overlap;

    vector < pair <uint, uint> > type_map_overlap;

    struct SampleInfo {

        double probability;
        const vector<bool> * type_p;
        string gamma;
    };

    void getLociTypes(const vector<bool> &, pair<vector<bool>, vector<bool> > *);

    template <typename T>
    struct ComparePair {
        
        bool operator() (T const& p1, T const& p2) const {

            return ((p1.first + p1.second) < (p2.first + p2.second));
        }
    };


  public:

    TypeContainer(mt19937*, AncestralTypeMaps &, pair<uint, uint> &);
    TypeContainer(mt19937*, pair<uint, uint> &);
    TypeContainer(const TypeContainer &);
    TypeContainer & operator=(const TypeContainer & in_copy);
    ~TypeContainer();

    bool isEmpty();
    uint getLineageCount();
    uint getLineageCountA();
    uint getLineageCountB();
    uint getLineageCountC();
    uint getLocusLengthA();
    uint getLocusLengthB();

    uint getTotalSequenceLength();
    uint getUniqueCount();

    uint getTypeOccurrenceC(Haplotype &);
    uint getTypeOccurrenceCa(Haplotype &);
    uint getTypeOccurrenceCb(Haplotype &);
    uint getTypeOccurrenceA(Haplotype &);
    uint getTypeOccurrenceB(Haplotype &);

    vector<Haplotype> getHaplotypesA();
    vector<Haplotype> getHaplotypesB();

    void addType(Haplotype &);
    void removeType(Haplotype &);
    void addTypeMap(AncestralTypeMaps &);
    vector<pair<uint, uint> > getTypeMapOverlap();
    pair<Haplotype, vector<double> >  sampleType(double, double, double);
    
    vector<vector<int> > getDifferences(Haplotype &);
    void calculateDifference(map < uint, uint > *, Haplotype &, TypeMap &);
};

#endif
