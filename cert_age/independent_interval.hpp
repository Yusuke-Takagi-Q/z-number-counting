// Independent dyadic interval arithmetic used for the certificate cross-checks.
#ifndef Z_CERTIFICATE_INDEPENDENT_INTERVAL_HPP
#define Z_CERTIFICATE_INDEPENDENT_INTERVAL_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace crosscheck {
using I = std::int64_t;
using Z = __int128_t;
constexpr I scale = I{1} << 50;
constexpr I lift = I{1} << 8;

inline I checked(Z x) {
    if (x < INT64_MIN || x > INT64_MAX) throw std::overflow_error("interval endpoint");
    return static_cast<I>(x);
}
inline I down(Z x, Z divisor) {
    return checked(x / divisor - (x % divisor < 0));
}
inline I up(Z x, Z divisor) {
    return checked(x / divisor + (x % divisor > 0));
}
struct Range { I lo = 0, hi = 0; };
struct Box { Range re, im; };

inline Range operator+(Range a, Range b) {
    return {checked(Z(a.lo) + b.lo), checked(Z(a.hi) + b.hi)};
}
inline Range operator-(Range a) {
    return {checked(-Z(a.hi)), checked(-Z(a.lo))};
}
inline Range multiply(Range a, Range b, Z divisor = scale) {
    const auto values = {Z(a.lo)*b.lo, Z(a.lo)*b.hi, Z(a.hi)*b.lo, Z(a.hi)*b.hi};
    return {down(std::min(values), divisor), up(std::max(values), divisor)};
}
inline Box operator+(Box a, Box b) {
    return {a.re + b.re, a.im + b.im};
}
inline Box operator*(Box a, Box b) {
    return {multiply(a.re, b.re) + -multiply(a.im, b.im),
            multiply(a.re, b.im) + multiply(a.im, b.re)};
}
inline Range widen_input(I lo, I hi) {
    return {checked(Z(lo)*lift), checked(Z(hi)*lift)};
}
inline Range shrink_output(Range a) {
    return {down(a.lo, lift), up(a.hi, lift)};
}
inline Box horner_exponential(Range angle) {
    // Degree 31, using 1+(it/1)(1+(it/2)(...)). The drivers
    // prove the real-argument remainder is at most 2^-42.
    const Range one{scale, scale};
    Box value{one, {}};
    for (I factorial_step = 31; factorial_step >= 1; --factorial_step) {
        const Z divisor = Z(scale)*factorial_step;
        value = {one + -multiply(angle, value.im, divisor),
                 multiply(angle, value.re, divisor)};
    }
    const Range error{-lift, lift};
    return {value.re+error, value.im+error};
}
inline I modulus_up(Box value) {
    const auto absolute = [](I i) { return i < 0 ? -Z(i) : Z(i); };
    const Z re = std::max(absolute(value.re.lo), absolute(value.re.hi));
    const Z im = std::max(absolute(value.im.lo), absolute(value.im.hi));
    if (re > (Z(1)<<60) || im > (Z(1)<<60)) throw std::overflow_error("modulus");
    const Z square = re*re + im*im;
    // The library square root is only a seed. Exact integer comparisons
    // determine the final result.
    const long double seed = std::sqrt(static_cast<long double>(square));
    if (!(seed >= 0 && seed < 0x1p62L)) throw std::overflow_error("square root seed");
    I root = static_cast<I>(seed);
    while (Z(root)*root < square) ++root;
    while (root > 0 && Z(root-1)*(root-1) >= square) --root;
    return root;
}
} // namespace crosscheck

#endif
