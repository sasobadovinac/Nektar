///////////////////////////////////////////////////////////////////////////////
//
// File: NektarExpMacros.h
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

#ifndef NEKTAR_LIBRARY_EXP_MACROS_H
#define NEKTAR_LIBRARY_EXP_MACROS_H

// The original macros used for registering before restructing the
// library. As the marcos can be difficult to set up and understand
// these are left in as an example. Going forward such macros could be
// used for other class registering. As such, it should probably live
// in a more general directory (and why it has generic name).

/*
#define NEKTAR_EXPL_DEC_INC1(C, R) \
std::string __register_##C##_##R            =
GetOperatorFactory().RegisterCreatorFunction( std::string( #C "_" #R "_Regular"
), &C##R::Create);

#define NEKTAR_EXPL_DEC_INC2(C, R) \
std::string __register_##C##_##R            =
GetOperatorFactory().RegisterCreatorFunction( std::string( #C "_" #R "_Regular"
), &C##R<>::Create); \ std::string __register_##C##_##R##_Deformed =
GetOperatorFactory().RegisterCreatorFunction( std::string( #C "_" #R "_Deformed"
), &C##R<true>::Create);
*/

#define NEKTAR_EXPL_DEC_INC_GEN1(C, R)                                         \
    std::string __register_##C##_##R =                                         \
        GetOperatorFactory().RegisterCreatorFunction(                          \
            std::string(#C "_" #R "_Regular"),                                 \
            &C##Template<LibUtilities::ShapeType::R>::Create);

#define NEKTAR_EXPL_DEC_INC_GEN2(C, R)                                         \
    std::string __register_##C##_##R =                                         \
        GetOperatorFactory().RegisterCreatorFunction(                          \
            std::string(#C "_" #R "_Regular"),                                 \
            &C##Template<LibUtilities::ShapeType::R>::Create);                 \
    std::string __register_##C##_##R##_Deformed =                              \
        GetOperatorFactory().RegisterCreatorFunction(                          \
            std::string(#C "_" #R "_Deformed"),                                \
            &C##Template<LibUtilities::ShapeType::R, true>::Create);

#endif
