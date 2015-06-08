
/*
Matrix3d.hpp - This file is part of the Coalescenator (v1.0.0)


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


#ifndef __coalescenator__Matrix3D_hpp
#define __coalescenator__Matrix3D_hpp

#include <vector>

#include "Utils.hpp"
#include "LikelihoodSamples.hpp"


using namespace std;


template <typename valueType>
class Matrix3d {

    private:

        const uint dim1_;
        const uint dim2_;
        const uint dim3_;

        vector<vector<vector<valueType> > > matrix_;

    public:

        Matrix3d(const uint, const uint, const uint, const uint);

        uint dim1();
        uint dim2();
        uint dim3();

        double getSummedElement(const uint i, const uint j, const uint k);
        valueType getElement(const uint i, const uint j, const uint k);

        void addValueToElement(const double, const uint, const uint, const uint);
        void addConstant(const double);
        void subtractConstant(const double);

        void addMatrixElementWise(Matrix3d<double> &);
        void addMatrixElementWise(Matrix3d<LikelihoodSamplesType> &);
        void logAddMatrixElementWise(Matrix3d<double> &);

        void multiplyValueAndAddMatrixElementWise(Matrix3d<double> &, const double);
        void addAndAddMatrixElementWise(Matrix3d<double> &, const double);
        void subtractMatrixElementWise(Matrix3d<double> &);

        void addSummedMatrixElementWise(Matrix3d<LikelihoodSamplesType> &);
};

#endif
