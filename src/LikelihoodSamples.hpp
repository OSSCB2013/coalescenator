
/*
LikelihoodSamples.hpp - This file is part of the Coalescenator (v1.0.0)


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


#ifndef __coalescenator__LikelihoodSamples_hpp
#define __coalescenator__LikelihoodSamples_hpp

#include <vector>
#include <set>

#include "Utils.hpp"

using namespace std;

class LikelihoodSamples {

    protected:

        multiset<double> likelihood_samples_;
        const uint num_samples;

    public:

        LikelihoodSamples(const uint);
        void operator+= (const double &);
        void operator+= (const LikelihoodSamples &);

        virtual double sumSamples() = 0;       

};


class LikelihoodSamplesSum : public LikelihoodSamples {

    public:

        LikelihoodSamplesSum(const uint);
        double sumSamples();
};


class LikelihoodSamplesKahan : public LikelihoodSamples {

    public:

        LikelihoodSamplesKahan(const uint);
        double sumSamples();
};

typedef LikelihoodSamplesSum LikelihoodSamplesType;


#endif
