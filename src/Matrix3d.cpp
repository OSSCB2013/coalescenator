
/*
Matrix3d.cpp - This file is part of the Coalescenator (v1.0.0)


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


#include <vector>
#include <assert.h>
#include <math.h>

#include "Matrix3d.hpp"
#include "Utils.hpp"
#include "LikelihoodSamples.hpp"

using namespace std;

template <typename valueType>
Matrix3d<valueType>::Matrix3d(const uint dim1_in, const uint dim2_in, const uint dim3_in, const uint num_samples) : dim1_(dim1_in), dim2_(dim2_in), dim3_(dim3_in) {
    
    matrix_ = vector < vector < vector < valueType > > > (dim1_, vector < vector < valueType > >(dim2_, vector < valueType >(dim3_, valueType(num_samples))));
}


template <typename valueType>
uint Matrix3d<valueType>::dim1() {

    return dim1_;
}


template <typename valueType>
uint Matrix3d<valueType>::dim2() {
    
    return dim2_;
}


template <typename valueType>
uint Matrix3d<valueType>::dim3() {

    return dim3_;
}


template <typename valueType> 
valueType Matrix3d<valueType>::getElement(const uint i, const uint j, const uint k) {

    return matrix_.at(i).at(j).at(k);
}


template <> 
double Matrix3d<LikelihoodSamplesType>::getSummedElement(const uint i, const uint j, const uint k) {

    return matrix_.at(i).at(j).at(k).sumSamples();
}


template <typename valueType> 
void Matrix3d<valueType>::addValueToElement(const double value, const uint i, const uint j, const uint k) {
    
    matrix_.at(i).at(j).at(k) += value; 
}


template <typename valueType> 
void Matrix3d<valueType>::addConstant(const double constant) {

    for (auto it_i = matrix_.begin(); it_i != matrix_.end(); it_i++) {

        for (auto it_j = it_i->begin(); it_j != it_i->end(); it_j++) {

            for (auto it_k = it_j->begin(); it_k != it_j->end(); it_k++) {
                
                *it_k += constant;
            }
        }
    }
}


template <> 
void Matrix3d<double>::subtractConstant(const double constant) {

    for (auto it_i = matrix_.begin(); it_i != matrix_.end(); it_i++) {

        for (auto it_j = it_i->begin(); it_j != it_i->end(); it_j++) {

            for (auto it_k = it_j->begin(); it_k != it_j->end(); it_k++) {
                
                *it_k -= constant;
            }
        }
    }
}


template <>
void Matrix3d<double>::addMatrixElementWise(Matrix3d<double> & matrix_in) {

    assert(dim1_ == matrix_in.dim1());
    assert(dim2_ == matrix_in.dim2());
    assert(dim3_ == matrix_in.dim3());
    
    for (auto i = 0; i < dim1_; i++) {

        for (auto j = 0; j < dim2_; j++) {

            for (auto k = 0; k < dim3_; k++) {
                
                matrix_.at(i).at(j).at(k) += matrix_in.getElement(i, j ,k);
            }
        }
    }
}


template <>
void Matrix3d<LikelihoodSamplesType>::addMatrixElementWise(Matrix3d<double> & matrix_in) {

    assert(dim1_ == matrix_in.dim1());
    assert(dim2_ == matrix_in.dim2());
    assert(dim3_ == matrix_in.dim3());
    
    for (auto i = 0; i < dim1_; i++) {

        for (auto j = 0; j < dim2_; j++) {

            for (auto k = 0; k < dim3_; k++) {
                
                matrix_.at(i).at(j).at(k) += matrix_in.getElement(i, j ,k);
            }
        }
    }
}


template <>
void Matrix3d<LikelihoodSamplesType>::addMatrixElementWise(Matrix3d<LikelihoodSamplesType> & matrix_in) {

    assert(dim1_ == matrix_in.dim1());
    assert(dim2_ == matrix_in.dim2());
    assert(dim3_ == matrix_in.dim3());
    
    for (auto i = 0; i < dim1_; i++) {

        for (auto j = 0; j < dim2_; j++) {

            for (auto k = 0; k < dim3_; k++) {
                
                matrix_.at(i).at(j).at(k) += matrix_in.getElement(i, j ,k);
            }
        }
    }
}


template <>
void Matrix3d<double>::logAddMatrixElementWise(Matrix3d<double> & matrix_in) {

    assert(dim1_ == matrix_in.dim1());
    assert(dim2_ == matrix_in.dim2());
    assert(dim3_ == matrix_in.dim3());
    
    for (auto i = 0; i < dim1_; i++) {

        for (auto j = 0; j < dim2_; j++) {

            for (auto k = 0; k < dim3_; k++) {
                
                matrix_.at(i).at(j).at(k) = logAddition(matrix_.at(i).at(j).at(k), matrix_in.getElement(i, j ,k));
            }
        }
    }
}



template <typename valueType>
void Matrix3d<valueType>::multiplyValueAndAddMatrixElementWise(Matrix3d<double> & matrix_in, const double value) {

    assert(dim1_ == matrix_in.dim1());
    assert(dim2_ == matrix_in.dim2());
    assert(dim3_ == matrix_in.dim3());
    
    for (auto i = 0; i < dim1_; i++) {

        for (auto j = 0; j < dim2_; j++) {

            for (auto k = 0; k < dim3_; k++) {
                
                matrix_.at(i).at(j).at(k) += matrix_in.getElement(i, j ,k) * value;
            }
        }
    }
}



template <typename valueType>
void Matrix3d<valueType>::addAndAddMatrixElementWise(Matrix3d<double> & matrix_in, const double value) {

    assert(dim1_ == matrix_in.dim1());
    assert(dim2_ == matrix_in.dim2());
    assert(dim3_ == matrix_in.dim3());
    
    for (auto i = 0; i < dim1_; i++) {

        for (auto j = 0; j < dim2_; j++) {

            for (auto k = 0; k < dim3_; k++) {
                
                matrix_.at(i).at(j).at(k) += matrix_in.getElement(i, j ,k) + value;
            }
        }
    }
}


template <>
void Matrix3d<double>::subtractMatrixElementWise(Matrix3d<double> & matrix_in) {

    assert(dim1_ == matrix_in.dim1());
    assert(dim2_ == matrix_in.dim2());
    assert(dim3_ == matrix_in.dim3());
    
    for (auto i = 0; i < dim1_; i++) {

        for (auto j = 0; j < dim2_; j++) {

            for (auto k = 0; k < dim3_; k++) {
                
                matrix_.at(i).at(j).at(k) -= matrix_in.getElement(i, j ,k);
            }
        }
    }
}

template <>
void Matrix3d<double>::addSummedMatrixElementWise(Matrix3d<LikelihoodSamplesType> & matrix_in) {
    
    assert(dim1_ == matrix_in.dim1());
    assert(dim2_ == matrix_in.dim2());
    assert(dim3_ == matrix_in.dim3());
    
    for (auto i = 0; i < dim1_; i++) {

        for (auto j = 0; j < dim2_; j++) {

            for (auto k = 0; k < dim3_; k++) {
                
                assert(matrix_.at(i).at(j).at(k) == 0);
                matrix_.at(i).at(j).at(k) += matrix_in.getSummedElement(i, j ,k);
            }
        }
    }   
}

template class Matrix3d<double>;
template class Matrix3d<LikelihoodSamplesType>;


