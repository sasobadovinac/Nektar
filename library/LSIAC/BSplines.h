///////////////////////////////////////////////////////////////////////////////
//
// File BSplines.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: BSplines definition
//
///////////////////////////////////////////////////////////////////////////////
#pragma once
#include "LSIACPostProcessor.h"
#include <iostream>
#include <vector>
using namespace std;

/// This class is the base class for all the B-Splines.
/** This class is useful when dynmaically creating object of its subclasses.
        All the BSplines need by the fitler would be subclasses of this filter.
*/
namespace Nektar
{
namespace LSIAC
{
class BSplines : public LSIACPostProcessor
{
protected:
public:
};
/*
class GeneralBSplines: public BSplines{
    // data
        public:
                vector<double>* m_knotVector;
                int m_Order;

        // functions
        public:
        GeneralBSplines(const vector<double> knots,const int order)
        {
                m_knotVector = new vector<double>(0);
                SetKnotVector(knots);
                SetOrder(order);
        }

        bool SetKnotVector(vector<double> knots)
        {
                this->m_knotVector->clear();
                this->m_knotVector->resize(knots.size());
                this->m_knotVector->assign(knots.begin(),knots.end());
                return true;
        }

        bool GetKnotVector ( vector<double> &knots)
        {
                knots = *m_knotVector;
                return true; }

        bool SetOrder(int order)
        {
                m_Order = order;
                return true; }

        int  GetOrder() const
        {return m_Order; }

        bool EvaluateBSplines(const  vector<double> t_pos,const vector<double>
knots, vector<double> &t_values) {return true; }

        bool EvaluateBSplines (const vector<double> t_pos, vector<double>
&t_values)const {return true; }

};

int main()
{
        cout << " Enterend into main functions " << endl;
        double mydoubles[] = { 0.0, 1.0, 2.0, 3.0, 4.0};
        vector<double> knots(mydoubles, mydoubles +
sizeof(mydoubles)/sizeof(double)); GeneralBSplines gsp(knots, 3); cout <<
gsp.m_knotVector->at(3)<< endl; cout <<  knots.size()<< endl; vector<double>
*knotsout; gsp.GetKnotVector(*knotsout); knotsout->at(3) = 10; cout <<
knotsout->size()<< endl; cout << knotsout->at(0)<< endl; cout <<
knotsout->at(1)<< endl; cout << knotsout->at(2)<< endl; cout <<
knotsout->at(3)<< endl; cout << gsp.m_knotVector->at(3)<< endl;

        return 1;
}
*/
} // namespace LSIAC
} // namespace Nektar
