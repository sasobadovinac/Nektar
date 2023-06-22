///////////////////////////////////////////////////////////////////////////////
//
// File: SwitchNodesPoints.h
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
// Description:
//
///////////////////////////////////////////////////////////////////////////////

// This code is the mombo switch statement that is used in mutliple
// operators. It uses preprocessor directives based on the shape
// type and dimension to limit the code inclusion.

// This header is included in the OperatorType.h file. Ideally the
// operator() method would be a template in the Operator base class
// but becasue the operator{1,23}D is both a function and template that
// need to be in the inherited class it is not possible.

#if defined(SHAPE_DIMENSION_1D)

{
    const int nm0 = m_basis[0]->GetNumModes();
    const int nq0 = m_basis[0]->GetNumPoints();

#if defined(SHAPE_TYPE_SEG)

    switch (nm0)
    {
        case 2:
            switch (nq0)
            {
                case 2:
                    operator1D<2, 2>(input, output);
                    break;
                case 3:
                    operator1D<2, 3>(input, output);
                    break;
                case 4:
                    operator1D<2, 4>(input, output);
                    break;
                default:
                    operator1D(input, output);
                    break;
            }
            break;
        case 3:
            switch (nq0)
            {
                case 3:
                    operator1D<3, 3>(input, output);
                    break;
                case 4:
                    operator1D<3, 4>(input, output);
                    break;
                case 5:
                    operator1D<3, 5>(input, output);
                    break;
                case 6:
                    operator1D<3, 6>(input, output);
                    break;
                default:
                    operator1D(input, output);
                    break;
            }
            break;
        case 4:
            switch (nq0)
            {
                case 4:
                    operator1D<4, 4>(input, output);
                    break;
                case 5:
                    operator1D<4, 5>(input, output);
                    break;
                case 6:
                    operator1D<4, 6>(input, output);
                    break;
                case 7:
                    operator1D<4, 7>(input, output);
                    break;
                case 8:
                    operator1D<4, 8>(input, output);
                    break;
                default:
                    operator1D(input, output);
                    break;
            }
            break;
        case 5:
            switch (nq0)
            {
                case 5:
                    operator1D<5, 5>(input, output);
                    break;
                case 6:
                    operator1D<5, 6>(input, output);
                    break;
                case 7:
                    operator1D<5, 7>(input, output);
                    break;
                case 8:
                    operator1D<5, 8>(input, output);
                    break;
                case 9:
                    operator1D<5, 9>(input, output);
                    break;
                case 10:
                    operator1D<5, 10>(input, output);
                    break;
                default:
                    operator1D(input, output);
                    break;
            }
            break;
        case 6:
            switch (nq0)
            {
                case 6:
                    operator1D<6, 6>(input, output);
                    break;
                case 7:
                    operator1D<6, 7>(input, output);
                    break;
                case 8:
                    operator1D<6, 8>(input, output);
                    break;
                case 9:
                    operator1D<6, 9>(input, output);
                    break;
                case 10:
                    operator1D<6, 10>(input, output);
                    break;
                case 11:
                    operator1D<6, 11>(input, output);
                    break;
                case 12:
                    operator1D<6, 12>(input, output);
                    break;
                default:
                    operator1D(input, output);
                    break;
            }
            break;
        case 7:
            switch (nq0)
            {
                case 7:
                    operator1D<7, 7>(input, output);
                    break;
                case 8:
                    operator1D<7, 8>(input, output);
                    break;
                case 9:
                    operator1D<7, 9>(input, output);
                    break;
                case 10:
                    operator1D<7, 10>(input, output);
                    break;
                case 11:
                    operator1D<7, 11>(input, output);
                    break;
                case 12:
                    operator1D<7, 12>(input, output);
                    break;
                case 13:
                    operator1D<7, 13>(input, output);
                    break;
                case 14:
                    operator1D<7, 14>(input, output);
                    break;
                default:
                    operator1D(input, output);
                    break;
            }
            break;
        case 8:
            switch (nq0)
            {
                case 8:
                    operator1D<8, 8>(input, output);
                    break;
                case 9:
                    operator1D<8, 9>(input, output);
                    break;
                case 10:
                    operator1D<8, 10>(input, output);
                    break;
                case 11:
                    operator1D<8, 11>(input, output);
                    break;
                case 12:
                    operator1D<8, 12>(input, output);
                    break;
                case 13:
                    operator1D<8, 13>(input, output);
                    break;
                case 14:
                    operator1D<8, 14>(input, output);
                    break;
                case 15:
                    operator1D<8, 15>(input, output);
                    break;
                case 16:
                    operator1D<8, 16>(input, output);
                    break;
                default:
                    operator1D(input, output);
                    break;
            }
            break;
            ;

        default:
            operator1D(input, output);
            break;
    }

#endif
}

#elif defined(SHAPE_DIMENSION_2D)

{
    const int nm0 = m_basis[0]->GetNumModes();
    const int nm1 = m_basis[1]->GetNumModes();

    const int nq0 = m_basis[0]->GetNumPoints();
    const int nq1 = m_basis[1]->GetNumPoints();

#if defined(SHAPE_TYPE_TRI)

    if (nm0 == nm1 && nq0 == nq1 + 1)
    {
        switch (nm0)
        {
            case 2:
                switch (nq0)
                {
                    case 3:
                        operator2D<2, 2, 3, 2>(input, output);
                        break;
                    case 4:
                        operator2D<2, 2, 4, 3>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 3:
                switch (nq0)
                {
                    case 4:
                        operator2D<3, 3, 4, 3>(input, output);
                        break;
                    case 5:
                        operator2D<3, 3, 5, 4>(input, output);
                        break;
                    case 6:
                        operator2D<3, 3, 6, 5>(input, output);
                        break;
                    case 7:
                        operator2D<3, 3, 7, 6>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 4:
                switch (nq0)
                {
                    case 5:
                        operator2D<4, 4, 5, 4>(input, output);
                        break;
                    case 6:
                        operator2D<4, 4, 6, 5>(input, output);
                        break;
                    case 7:
                        operator2D<4, 4, 7, 6>(input, output);
                        break;
                    case 8:
                        operator2D<4, 4, 8, 7>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 5:
                switch (nq0)
                {
                    case 6:
                        operator2D<5, 5, 6, 5>(input, output);
                        break;
                    case 7:
                        operator2D<5, 5, 7, 6>(input, output);
                        break;
                    case 8:
                        operator2D<5, 5, 8, 7>(input, output);
                        break;
                    case 9:
                        operator2D<5, 5, 9, 8>(input, output);
                        break;
                    case 10:
                        operator2D<5, 5, 10, 9>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 6:
                switch (nq0)
                {
                    case 7:
                        operator2D<6, 6, 7, 6>(input, output);
                        break;
                    case 8:
                        operator2D<6, 6, 8, 7>(input, output);
                        break;
                    case 9:
                        operator2D<6, 6, 9, 8>(input, output);
                        break;
                    case 10:
                        operator2D<6, 6, 10, 9>(input, output);
                        break;
                    case 11:
                        operator2D<6, 6, 11, 10>(input, output);
                        break;
                    case 12:
                        operator2D<6, 6, 12, 11>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 7:
                switch (nq0)
                {
                    case 8:
                        operator2D<7, 7, 8, 7>(input, output);
                        break;
                    case 9:
                        operator2D<7, 7, 9, 8>(input, output);
                        break;
                    case 10:
                        operator2D<7, 7, 10, 9>(input, output);
                        break;
                    case 11:
                        operator2D<7, 7, 11, 10>(input, output);
                        break;
                    case 12:
                        operator2D<7, 7, 12, 11>(input, output);
                        break;
                    case 13:
                        operator2D<7, 7, 13, 12>(input, output);
                        break;
                    case 14:
                        operator2D<7, 7, 14, 13>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 8:
                switch (nq0)
                {
                    case 9:
                        operator2D<8, 8, 9, 8>(input, output);
                        break;
                    case 10:
                        operator2D<8, 8, 10, 9>(input, output);
                        break;
                    case 11:
                        operator2D<8, 8, 11, 10>(input, output);
                        break;
                    case 12:
                        operator2D<8, 8, 12, 11>(input, output);
                        break;
                    case 13:
                        operator2D<8, 8, 13, 12>(input, output);
                        break;
                    case 14:
                        operator2D<8, 8, 14, 13>(input, output);
                        break;
                    case 15:
                        operator2D<8, 8, 15, 14>(input, output);
                        break;
                    case 16:
                        operator2D<8, 8, 16, 15>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            default:
                operator2D(input, output);
                break;
        }
    }

#elif defined(SHAPE_TYPE_QUAD)

    if (nm0 == nm1 && nq0 == nq1)
    {
        switch (nm0)
        {
            case 2:
                switch (nq0)
                {
                    case 2:
                        operator2D<2, 2, 2, 2>(input, output);
                        break;
                    case 3:
                        operator2D<2, 2, 3, 3>(input, output);
                        break;
                    case 4:
                        operator2D<2, 2, 4, 4>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 3:
                switch (nq0)
                {
                    case 3:
                        operator2D<3, 3, 3, 3>(input, output);
                        break;
                    case 4:
                        operator2D<3, 3, 4, 4>(input, output);
                        break;
                    case 5:
                        operator2D<3, 3, 5, 5>(input, output);
                        break;
                    case 6:
                        operator2D<3, 3, 6, 6>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 4:
                switch (nq0)
                {
                    case 4:
                        operator2D<4, 4, 4, 4>(input, output);
                        break;
                    case 5:
                        operator2D<4, 4, 5, 5>(input, output);
                        break;
                    case 6:
                        operator2D<4, 4, 6, 6>(input, output);
                        break;
                    case 7:
                        operator2D<4, 4, 7, 7>(input, output);
                        break;
                    case 8:
                        operator2D<4, 4, 8, 8>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 5:
                switch (nq0)
                {
                    case 5:
                        operator2D<5, 5, 5, 5>(input, output);
                        break;
                    case 6:
                        operator2D<5, 5, 6, 6>(input, output);
                        break;
                    case 7:
                        operator2D<5, 5, 7, 7>(input, output);
                        break;
                    case 8:
                        operator2D<5, 5, 8, 8>(input, output);
                        break;
                    case 9:
                        operator2D<5, 5, 9, 9>(input, output);
                        break;
                    case 10:
                        operator2D<5, 5, 10, 10>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 6:
                switch (nq0)
                {
                    case 6:
                        operator2D<6, 6, 6, 6>(input, output);
                        break;
                    case 7:
                        operator2D<6, 6, 7, 7>(input, output);
                        break;
                    case 8:
                        operator2D<6, 6, 8, 8>(input, output);
                        break;
                    case 9:
                        operator2D<6, 6, 9, 9>(input, output);
                        break;
                    case 10:
                        operator2D<6, 6, 10, 10>(input, output);
                        break;
                    case 11:
                        operator2D<6, 6, 11, 11>(input, output);
                        break;
                    case 12:
                        operator2D<6, 6, 12, 12>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 7:
                switch (nq0)
                {
                    case 7:
                        operator2D<7, 7, 7, 7>(input, output);
                        break;
                    case 8:
                        operator2D<7, 7, 8, 8>(input, output);
                        break;
                    case 9:
                        operator2D<7, 7, 9, 9>(input, output);
                        break;
                    case 10:
                        operator2D<7, 7, 10, 10>(input, output);
                        break;
                    case 11:
                        operator2D<7, 7, 11, 11>(input, output);
                        break;
                    case 12:
                        operator2D<7, 7, 12, 12>(input, output);
                        break;
                    case 13:
                        operator2D<7, 7, 13, 13>(input, output);
                        break;
                    case 14:
                        operator2D<7, 7, 14, 14>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            case 8:
                switch (nq0)
                {
                    case 8:
                        operator2D<8, 8, 8, 8>(input, output);
                        break;
                    case 9:
                        operator2D<8, 8, 9, 9>(input, output);
                        break;
                    case 10:
                        operator2D<8, 8, 10, 10>(input, output);
                        break;
                    case 11:
                        operator2D<8, 8, 11, 11>(input, output);
                        break;
                    case 12:
                        operator2D<8, 8, 12, 12>(input, output);
                        break;
                    case 13:
                        operator2D<8, 8, 13, 13>(input, output);
                        break;
                    case 14:
                        operator2D<8, 8, 14, 14>(input, output);
                        break;
                    case 15:
                        operator2D<8, 8, 15, 15>(input, output);
                        break;
                    case 16:
                        operator2D<8, 8, 16, 16>(input, output);
                        break;
                    default:
                        operator2D(input, output);
                        break;
                }
                break;
            default:
                operator2D(input, output);
                break;
        }
    }

#endif // SHAPE_TYPE

    else
    {
        operator2D(input, output);
    }
}

#elif defined(SHAPE_DIMENSION_3D)

{
    const int nm0 = m_basis[0]->GetNumModes();
    const int nm1 = m_basis[1]->GetNumModes();
    const int nm2 = m_basis[2]->GetNumModes();

    const int nq0 = m_basis[0]->GetNumPoints();
    const int nq1 = m_basis[1]->GetNumPoints();
    const int nq2 = m_basis[2]->GetNumPoints();

#if defined(SHAPE_TYPE_HEX)

    if (nm0 == nm1 && nm0 == nm2 && nq0 == nq1 && nq0 == nq2)
    {
        switch (nm0)
        {
            case 2:
                switch (nq0)
                {
                    case 2:
                        operator3D<2, 2, 2, 2, 2, 2>(input, output);
                        break;
                    case 3:
                        operator3D<2, 2, 2, 3, 3, 3>(input, output);
                        break;
                    case 4:
                        operator3D<2, 2, 2, 4, 4, 4>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 3:
                switch (nq0)
                {
                    case 3:
                        operator3D<3, 3, 3, 3, 3, 3>(input, output);
                        break;
                    case 4:
                        operator3D<3, 3, 3, 4, 4, 4>(input, output);
                        break;
                    case 5:
                        operator3D<3, 3, 3, 5, 5, 5>(input, output);
                        break;
                    case 6:
                        operator3D<3, 3, 3, 6, 6, 6>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 4:
                switch (nq0)
                {
                    case 4:
                        operator3D<4, 4, 4, 4, 4, 4>(input, output);
                        break;
                    case 5:
                        operator3D<4, 4, 4, 5, 5, 5>(input, output);
                        break;
                    case 6:
                        operator3D<4, 4, 4, 6, 6, 6>(input, output);
                        break;
                    case 7:
                        operator3D<4, 4, 4, 7, 7, 7>(input, output);
                        break;
                    case 8:
                        operator3D<4, 4, 4, 8, 8, 8>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 5:
                switch (nq0)
                {
                    case 5:
                        operator3D<5, 5, 5, 5, 5, 5>(input, output);
                        break;
                    case 6:
                        operator3D<5, 5, 5, 6, 6, 6>(input, output);
                        break;
                    case 7:
                        operator3D<5, 5, 5, 7, 7, 7>(input, output);
                        break;
                    case 8:
                        operator3D<5, 5, 5, 8, 8, 8>(input, output);
                        break;
                    case 9:
                        operator3D<5, 5, 5, 9, 9, 9>(input, output);
                        break;
                    case 10:
                        operator3D<5, 5, 5, 10, 10, 10>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 6:
                switch (nq0)
                {
                    case 6:
                        operator3D<6, 6, 6, 6, 6, 6>(input, output);
                        break;
                    case 7:
                        operator3D<6, 6, 6, 7, 7, 7>(input, output);
                        break;
                    case 8:
                        operator3D<6, 6, 6, 8, 8, 8>(input, output);
                        break;
                    case 9:
                        operator3D<6, 6, 6, 9, 9, 9>(input, output);
                        break;
                    case 10:
                        operator3D<6, 6, 6, 10, 10, 10>(input, output);
                        break;
                    case 11:
                        operator3D<6, 6, 6, 11, 11, 11>(input, output);
                        break;
                    case 12:
                        operator3D<6, 6, 6, 12, 12, 12>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 7:
                switch (nq0)
                {
                    case 7:
                        operator3D<7, 7, 7, 7, 7, 7>(input, output);
                        break;
                    case 8:
                        operator3D<7, 7, 7, 8, 8, 8>(input, output);
                        break;
                    case 9:
                        operator3D<7, 7, 7, 9, 9, 9>(input, output);
                        break;
                    case 10:
                        operator3D<7, 7, 7, 10, 10, 10>(input, output);
                        break;
                    case 11:
                        operator3D<7, 7, 7, 11, 11, 11>(input, output);
                        break;
                    case 12:
                        operator3D<7, 7, 7, 12, 12, 12>(input, output);
                        break;
                    case 13:
                        operator3D<7, 7, 7, 13, 13, 13>(input, output);
                        break;
                    case 14:
                        operator3D<7, 7, 7, 14, 14, 14>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 8:
                switch (nq0)
                {
                    case 8:
                        operator3D<8, 8, 8, 8, 8, 8>(input, output);
                        break;
                    case 9:
                        operator3D<8, 8, 8, 9, 9, 9>(input, output);
                        break;
                    case 10:
                        operator3D<8, 8, 8, 10, 10, 10>(input, output);
                        break;
                    case 11:
                        operator3D<8, 8, 8, 11, 11, 11>(input, output);
                        break;
                    case 12:
                        operator3D<8, 8, 8, 12, 12, 12>(input, output);
                        break;
                    case 13:
                        operator3D<8, 8, 8, 13, 13, 13>(input, output);
                        break;
                    case 14:
                        operator3D<8, 8, 8, 14, 14, 14>(input, output);
                        break;
                    case 15:
                        operator3D<8, 8, 8, 15, 15, 15>(input, output);
                        break;
                    case 16:
                        operator3D<8, 8, 8, 16, 16, 16>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            default:
                operator3D(input, output);
                break;
        }
    }

#elif defined(SHAPE_TYPE_TET)

    if (nm0 == nm1 && nm0 == nm2 && nq0 == nq1 + 1 && nq0 == nq2 + 1)
    {
        switch (nm0)
        {
            case 2:
                switch (nq0)
                {
                    case 3:
                        operator3D<2, 2, 2, 3, 2, 2>(input, output);
                        break;
                    case 4:
                        operator3D<2, 2, 2, 4, 3, 3>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 3:
                switch (nq0)
                {
                    case 4:
                        operator3D<3, 3, 3, 4, 3, 3>(input, output);
                        break;
                    case 5:
                        operator3D<3, 3, 3, 5, 4, 4>(input, output);
                        break;
                    case 6:
                        operator3D<3, 3, 3, 6, 5, 5>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 4:
                switch (nq0)
                {
                    case 5:
                        operator3D<4, 4, 4, 5, 4, 4>(input, output);
                        break;
                    case 6:
                        operator3D<4, 4, 4, 6, 5, 5>(input, output);
                        break;
                    case 7:
                        operator3D<4, 4, 4, 7, 6, 6>(input, output);
                        break;
                    case 8:
                        operator3D<4, 4, 4, 8, 7, 7>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 5:
                switch (nq0)
                {
                    case 6:
                        operator3D<5, 5, 5, 6, 5, 5>(input, output);
                        break;
                    case 7:
                        operator3D<5, 5, 5, 7, 6, 6>(input, output);
                        break;
                    case 8:
                        operator3D<5, 5, 5, 8, 7, 7>(input, output);
                        break;
                    case 9:
                        operator3D<5, 5, 5, 9, 8, 8>(input, output);
                        break;
                    case 10:
                        operator3D<5, 5, 5, 10, 9, 9>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 6:
                switch (nq0)
                {
                    case 7:
                        operator3D<6, 6, 6, 7, 6, 6>(input, output);
                        break;
                    case 8:
                        operator3D<6, 6, 6, 8, 7, 7>(input, output);
                        break;
                    case 9:
                        operator3D<6, 6, 6, 9, 8, 8>(input, output);
                        break;
                    case 10:
                        operator3D<6, 6, 6, 10, 9, 9>(input, output);
                        break;
                    case 11:
                        operator3D<6, 6, 6, 11, 10, 10>(input, output);
                        break;
                    case 12:
                        operator3D<6, 6, 6, 12, 11, 11>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 7:
                switch (nq0)
                {
                    case 8:
                        operator3D<7, 7, 7, 8, 7, 7>(input, output);
                        break;
                    case 9:
                        operator3D<7, 7, 7, 9, 8, 8>(input, output);
                        break;
                    case 10:
                        operator3D<7, 7, 7, 10, 9, 9>(input, output);
                        break;
                    case 11:
                        operator3D<7, 7, 7, 11, 10, 10>(input, output);
                        break;
                    case 12:
                        operator3D<7, 7, 7, 12, 11, 11>(input, output);
                        break;
                    case 13:
                        operator3D<7, 7, 7, 13, 12, 12>(input, output);
                        break;
                    case 14:
                        operator3D<7, 7, 7, 14, 13, 13>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 8:
                switch (nq0)
                {
                    case 9:
                        operator3D<8, 8, 8, 9, 8, 8>(input, output);
                        break;
                    case 10:
                        operator3D<8, 8, 8, 10, 9, 9>(input, output);
                        break;
                    case 11:
                        operator3D<8, 8, 8, 11, 10, 10>(input, output);
                        break;
                    case 12:
                        operator3D<8, 8, 8, 12, 11, 11>(input, output);
                        break;
                    case 13:
                        operator3D<8, 8, 8, 13, 12, 12>(input, output);
                        break;
                    case 14:
                        operator3D<8, 8, 8, 14, 13, 13>(input, output);
                        break;
                    case 15:
                        operator3D<8, 8, 8, 15, 14, 14>(input, output);
                        break;
                    case 16:
                        operator3D<8, 8, 8, 16, 15, 15>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            default:
                operator3D(input, output);
                break;
        }
    }

#elif defined(SHAPE_TYPE_PYR) || defined(SHAPE_TYPE_PRISM)

    if (nm0 == nm1 && nm0 == nm2 && nq0 == nq1 && nq0 == nq2 + 1)
    {
        switch (nm0)
        {
            case 2:
                switch (nq0)
                {
                    case 3:
                        operator3D<2, 2, 2, 3, 3, 2>(input, output);
                        break;
                    case 4:
                        operator3D<2, 2, 2, 4, 4, 3>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 3:
                switch (nq0)
                {
                    case 4:
                        operator3D<3, 3, 3, 4, 4, 3>(input, output);
                        break;
                    case 5:
                        operator3D<3, 3, 3, 5, 5, 4>(input, output);
                        break;
                    case 6:
                        operator3D<3, 3, 3, 6, 6, 5>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 4:
                switch (nq0)
                {
                    case 5:
                        operator3D<4, 4, 4, 5, 5, 4>(input, output);
                        break;
                    case 6:
                        operator3D<4, 4, 4, 6, 6, 5>(input, output);
                        break;
                    case 7:
                        operator3D<4, 4, 4, 7, 7, 6>(input, output);
                        break;
                    case 8:
                        operator3D<4, 4, 4, 8, 8, 7>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 5:
                switch (nq0)
                {
                    case 6:
                        operator3D<5, 5, 5, 6, 6, 5>(input, output);
                        break;
                    case 7:
                        operator3D<5, 5, 5, 7, 7, 6>(input, output);
                        break;
                    case 8:
                        operator3D<5, 5, 5, 8, 8, 7>(input, output);
                        break;
                    case 9:
                        operator3D<5, 5, 5, 9, 9, 8>(input, output);
                        break;
                    case 10:
                        operator3D<5, 5, 5, 10, 10, 9>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 6:
                switch (nq0)
                {
                    case 7:
                        operator3D<6, 6, 6, 7, 7, 6>(input, output);
                        break;
                    case 8:
                        operator3D<6, 6, 6, 8, 8, 7>(input, output);
                        break;
                    case 9:
                        operator3D<6, 6, 6, 9, 9, 8>(input, output);
                        break;
                    case 10:
                        operator3D<6, 6, 6, 10, 10, 9>(input, output);
                        break;
                    case 11:
                        operator3D<6, 6, 6, 11, 11, 10>(input, output);
                        break;
                    case 12:
                        operator3D<6, 6, 6, 12, 12, 11>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 7:
                switch (nq0)
                {
                    case 8:
                        operator3D<7, 7, 7, 8, 8, 7>(input, output);
                        break;
                    case 9:
                        operator3D<7, 7, 7, 9, 9, 8>(input, output);
                        break;
                    case 10:
                        operator3D<7, 7, 7, 10, 10, 9>(input, output);
                        break;
                    case 11:
                        operator3D<7, 7, 7, 11, 11, 10>(input, output);
                        break;
                    case 12:
                        operator3D<7, 7, 7, 12, 12, 11>(input, output);
                        break;
                    case 13:
                        operator3D<7, 7, 7, 13, 13, 12>(input, output);
                        break;
                    case 14:
                        operator3D<7, 7, 7, 14, 14, 13>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            case 8:
                switch (nq0)
                {
                    case 9:
                        operator3D<8, 8, 8, 9, 9, 8>(input, output);
                        break;
                    case 10:
                        operator3D<8, 8, 8, 10, 10, 9>(input, output);
                        break;
                    case 11:
                        operator3D<8, 8, 8, 11, 11, 10>(input, output);
                        break;
                    case 12:
                        operator3D<8, 8, 8, 12, 12, 11>(input, output);
                        break;
                    case 13:
                        operator3D<8, 8, 8, 13, 13, 12>(input, output);
                        break;
                    case 14:
                        operator3D<8, 8, 8, 14, 14, 13>(input, output);
                        break;
                    case 15:
                        operator3D<8, 8, 8, 15, 15, 14>(input, output);
                        break;
                    case 16:
                        operator3D<8, 8, 8, 16, 16, 15>(input, output);
                        break;
                    default:
                        operator3D(input, output);
                        break;
                }
                break;
            default:
                operator3D(input, output);
                break;
        }
    }

#endif // SHAPE_TYPE

    else
    {
        operator3D(input, output);
    }
}

#endif // SHAPE_DIMENSION
