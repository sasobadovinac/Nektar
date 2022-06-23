// This code is the mombo switch statement that is used in mutliple
// operators. It uses preprocessor directives based on the shape
// type and dimension to limit the code inclusion.

// This header is included in the OperatorType.h file. Ideally the
// operator() method would be a template in the Operator base class
// but becasue the operator{1,23}D is both a function and template that
// need to be in the inherited class it is not possible.

#if defined(SHAPE_DIMENSION_1D)

    {
        const int nq0 = m_basis[0]->GetNumPoints();

#if defined(SHAPE_TYPE_SEG)

        switch(nq0)
        {
        case 2:
            operator1D<2>(input, output); break;
        case 3:
            operator1D<3>(input, output); break;
        case 4:
            operator1D<4>(input, output); break;
        case 5:
            operator1D<5>(input, output); break;
        case 6:
            operator1D<6>(input, output); break;
        case 7:
            operator1D<7>(input, output); break;
        case 8:
            operator1D<8>(input, output); break;
        case 9:
            operator1D<9>(input, output); break;
        case 10:
            operator1D<10>(input, output); break;
        default:
            operator1D(input, output); break;
        }

#endif

    }

#elif defined(SHAPE_DIMENSION_2D)

    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();


#if defined(SHAPE_TYPE_TRI)

        if(nq0 == nq1+1)
        {
            switch(m_basis[0]->GetNumPoints())
            {
            case 2:
                operator2D<2,1>(input, output); break;
            case 3:
                operator2D<3,2>(input, output); break;
            case 4:
                operator2D<4,3>(input, output); break;
            case 5:
                operator2D<5,4>(input, output); break;
            case 6:
                operator2D<6,5>(input, output); break;
            case 7:
                operator2D<7,6>(input, output); break;
            case 8:
                operator2D<8,7>(input, output); break;
            case 9:
                operator2D<9,8>(input, output); break;
            case 10:
                operator2D<10,9>(input, output); break;
            default:
                operator2D(input, output); break;
            }
        }

#elif defined(SHAPE_TYPE_QUAD)

        if(nq0 == nq1)
        {
            switch(m_basis[0]->GetNumPoints())
            {
            case 2:
                operator2D<2,2>(input, output); break;
            case 3:
                operator2D<3,3>(input, output); break;
            case 4:
                operator2D<4,4>(input, output); break;
            case 5:
                operator2D<5,5>(input, output); break;
            case 6:
                operator2D<6,6>(input, output); break;
            case 7:
                operator2D<7,7>(input, output); break;
            case 8:
                operator2D<8,8>(input, output); break;
            case 9:
                operator2D<9,9>(input, output); break;
            case 10:
                operator2D<10,10>(input, output); break;
            default:
                operator2D(input, output); break;
            }
        }

#endif  // SHAPE_TYPE

        else
        {
            operator2D(input, output);
        }
    }

#elif defined(SHAPE_DIMENSION_3D)

    {
        const int nq0 = m_basis[0]->GetNumPoints();
        const int nq1 = m_basis[1]->GetNumPoints();
        const int nq2 = m_basis[2]->GetNumPoints();

#if defined(SHAPE_TYPE_HEX)

        if(nq0 == nq1 && nq0 == nq2)
        {
            switch(m_basis[0]->GetNumPoints())
            {
            case 2:
                operator3D<2,2,2>(input, output); break;
            case 3:
                operator3D<3,3,3>(input, output); break;
            case 4:
                operator3D<4,4,4>(input, output); break;
            case 5:
                operator3D<5,5,5>(input, output); break;
            case 6:
                operator3D<6,6,6>(input, output); break;
            case 7:
                operator3D<7,7,7>(input, output); break;
            case 8:
                operator3D<8,8,8>(input, output); break;
            case 9:
                operator3D<9,9,9>(input, output); break;
            case 10:
                operator3D<10,10,10>(input, output); break;
            case 11:
                operator3D<11,11,11>(input, output); break;
            case 12:
                operator3D<12,12,12>(input, output); break;
            case 13:
                operator3D<13,13,13>(input, output); break;
            case 14:
                operator3D<14,14,14>(input, output); break;
            case 15:
                operator3D<15,15,15>(input, output); break;
            case 16:
                operator3D<16,16,16>(input, output); break;
            default:
                operator3D(input, output); break;
            }
        }

#elif defined(SHAPE_TYPE_TET)

        if(nq0 == nq1+1 && nq0 == nq2+1)
        {
            switch(nq0)
            {
            case 3:
                operator3D<3,2,2>(input, output); break;
            case 4:
                operator3D<4,3,3>(input, output); break;
            case 5:
                operator3D<5,4,4>(input, output); break;
            case 6:
                operator3D<6,5,5>(input, output); break;
            case 7:
                operator3D<7,6,6>(input, output); break;
            case 8:
                operator3D<8,7,7>(input, output); break;
            case 9:
                operator3D<9,8,8>(input, output); break;
            case 10:
                operator3D<10,9,9>(input, output); break;
            case 11:
                operator3D<11,10,10>(input, output); break;
            case 12:
                operator3D<12,11,11>(input, output); break;
            case 13:
                operator3D<13,12,12>(input, output); break;
            case 14:
                operator3D<14,13,13>(input, output); break;
            case 15:
                operator3D<15,14,14>(input, output); break;
            case 16:
                operator3D<16,15,15>(input, output); break;
            default:
                operator3D(input, output); break;
            }
        }

#elif defined(SHAPE_TYPE_PYR) || defined(SHAPE_TYPE_PRISM)

        if(nq0 == nq1 && nq0 == nq2+1)
        {
            switch(m_basis[0]->GetNumPoints())
            {
            case 2:
                operator3D<2,2,1>(input, output); break;
            case 3:
                operator3D<3,3,2>(input, output); break;
            case 4:
                operator3D<4,4,3>(input, output); break;
            case 5:
                operator3D<5,5,4>(input, output); break;
            case 6:
                operator3D<6,6,5>(input, output); break;
            case 7:
                operator3D<7,7,6>(input, output); break;
            case 8:
                operator3D<8,8,7>(input, output); break;
            case 9:
                operator3D<9,9,8>(input, output); break;
            case 10:
                operator3D<10,10,9>(input, output); break;
            case 11:
                operator3D<11,11,10>(input, output); break;
            case 12:
                operator3D<12,12,11>(input, output); break;
            case 13:
                operator3D<13,13,12>(input, output); break;
            case 14:
                operator3D<14,14,13>(input, output); break;
            case 15:
                operator3D<15,15,14>(input, output); break;
            case 16:
                operator3D<16,16,15>(input, output); break;
            default:
                operator3D(input, output); break;
            }
        }

#endif  // SHAPE_TYPE

        else
        {
            operator3D(input, output);
        }
    }

#endif  // SHAPE_DIMENSION
