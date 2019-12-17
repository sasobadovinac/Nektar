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
