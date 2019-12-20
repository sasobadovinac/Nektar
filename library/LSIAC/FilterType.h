///////////////////////////////////////////////////////////////////////////////
//
// File FilterType.h
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
// Description: FilterType definition
//
///////////////////////////////////////////////////////////////////////////////
#pragma once

namespace Nektar
{
namespace LSIAC
{
namespace SIACUtilities
{
//! This enum specifies the kind of filter to be used by the object.
/*!  The user can specify any type of smoothing or derivative filter.
 *  All the properties of memeber depends on this enum. This enum cannot be
 * changed after the object has been initialized. A new enum needs to be added
 * to the list to create any new type of filter.
 */
enum FilterType
{
    eSYM_2kp1_1SIDED_2kp2, /**< Default, Xli 2k+1 CentralBspl and 1 GeneralBSlp
                              only at borders appropriately*/
    eSYM_DER_2kp1_1SIDED_2kp1, /**< Simple Derivative filter*/
    eSYM_2kp1, /**< Will apply filter only when symmetric filter is possible. */
    eSYM_2kp1_1SIDED_2kp1, /**< Simplest form of SIAC and oneSided filter. */
    eSYM_2kp1_1SIDED_4kp1, /**< SRV's 4k+1 CentralBspl at boundaries.*/
    eSYM_DER_2kp1_1SIDED_2kp2, /**< SRV's Derivative filter */
    eSYM_DER_2kp1_1SIDED_4kp1, /**< Xli's Derivative fitler */
    eSYM_4kp1,
    eSYM_NDER_2kp1_1SIDED_2kp1, /**< This is for derivative Filter but not
                                   applying Divided Difference> */
    eSYM_UNEVEN_2kp1, /**< This is for new type of filter with uneven knots> */
    eNONE
};

struct SMOOTHIE_Status
{
    int count;
    int cancelled;
    int SOURCE;
    int TAG;
    int ERROR;
};
} // namespace SIACUtilities

} // namespace LSIAC
} // namespace Nektar
