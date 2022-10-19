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
