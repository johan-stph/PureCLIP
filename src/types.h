// ======================================================================
// PureCLIP — project-wide type definitions
// ======================================================================
// Float precision: double by default; long double with -DPURECLIP_HIGH_PRECISION
// On ARM64 (Apple Silicon), long double is software-emulated (~3× slower).
//
// NOTE: Switching to double shifts the sensitivity/specificity curve slightly
// (the GSL simplex converges to a marginally different local minimum).
// This fork produces ~63 crosslink sites on chr21 vs the original 86.
// To match the published output exactly, build with:
//   cmake -DPURECLIP_HIGH_PRECISION=ON -DPURECLIP_OMP_SCHEDULE=dynamic,1
// ======================================================================

#ifndef PURECLIP_TYPES_H_
#define PURECLIP_TYPES_H_

#include <cfloat>
#include <cmath>
#include <limits>
#include <type_traits>

#ifdef PURECLIP_HIGH_PRECISION
using Float = long double;
#else
using Float = double;
#endif

// Constexpr-friendly numeric limits
namespace pureclip {
inline constexpr Float float_min() noexcept {
    if constexpr (std::is_same_v<Float, long double>)
        return LDBL_MIN;
    else
        return DBL_MIN;
}
inline constexpr Float float_quiet_nan() noexcept {
    return std::numeric_limits<Float>::quiet_NaN();
}
inline constexpr int float_min_10_exp() noexcept {
    if constexpr (std::is_same_v<Float, long double>)
        return LDBL_MIN_10_EXP;
    else
        return DBL_MIN_10_EXP;
}
inline constexpr int float_digits10() noexcept {
    if constexpr (std::is_same_v<Float, long double>)
        return LDBL_DIG;
    else
        return DBL_DIG;
}
} // namespace pureclip

#endif // PURECLIP_TYPES_H_
