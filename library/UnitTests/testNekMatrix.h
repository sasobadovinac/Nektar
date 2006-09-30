///////////////////////////////////////////////////////////////////////////////
//
// File: testNekMatrix.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Tests NekMatrix functionality.
//
///////////////////////////////////////////////////////////////////////////////


#ifndef NEKTAR_UNIT_TESTS_TEST_NEK_MATRIX_H
#define NEKTAR_UNIT_TESTS_TEST_NEK_MATRIX_H

namespace Nektar
{
    namespace UnitTests
    {
        void testNekMatrixConstruction();
        void testNekMatrixAccess();
        void testNekMatrixBasicMath();
        void testNekMatrixFullDiagonalOperations();
        void testDiagonalMatrix();
        void testUserManagedMatrixData();
        void testBlockMatrices();
    }
}

#endif // NEKTAR_UNIT_TESTS_TEST_NEK_MATRIX_H


/**
    $Log: testNekMatrix.h,v $
    Revision 1.5  2006/08/25 01:38:59  bnelson
    no message

    Revision 1.4  2006/08/25 01:36:25  bnelson
    no message

    Revision 1.3  2006/08/14 02:35:45  bnelson
    Added many LinearAlgebra tests

    Revision 1.2  2006/05/31 04:19:37  bnelson
    Removed a test for invalid access to a matrix.

    Revision 1.1  2006/05/07 21:10:10  bnelson
    *** empty log message ***

**/

