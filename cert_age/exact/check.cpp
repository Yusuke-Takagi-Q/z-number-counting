// Fixed-point checker for the rational enclosures produced by prepare.py.
#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

using Int = std::int64_t;
using Wide = __int128_t;
constexpr Int scale = Int{1} << 42;
#ifndef CERT_WIDTH
#define CERT_WIDTH 7
#endif
constexpr std::size_t width = CERT_WIDTH;
static_assert(3 <= width && width <= 7);
constexpr unsigned masks = 1u << width;

Int narrow(Wide x) {
    assert(x >= std::numeric_limits<Int>::min() &&
           x <= std::numeric_limits<Int>::max());
    return static_cast<Int>(x);
}

Int floor_div(Wide x, Int divisor) {
    assert(divisor > 0);
    return narrow(x/divisor - (x%divisor < 0));
}

Int ceil_div(Wide x, Int divisor) {
    return -floor_div(-x, divisor);
}

struct Interval {
    Int lo = 0;
    Int hi = 0;
};

Interval operator+(Interval a, Interval b) {
    return {narrow(Wide(a.lo)+b.lo), narrow(Wide(a.hi)+b.hi)};
}

Interval operator-(Interval a) {
    return {narrow(-Wide(a.hi)), narrow(-Wide(a.lo))};
}

Interval operator*(Interval a, Interval b) {
    const auto products = {Wide(a.lo)*b.lo, Wide(a.lo)*b.hi,
                           Wide(a.hi)*b.lo, Wide(a.hi)*b.hi};
    return {floor_div(std::min(products), scale), ceil_div(std::max(products), scale)};
}

Int mul_up(Int a, Int b) {
    assert(a >= 0 && b >= 0);
    return ceil_div(Wide(a)*b, scale);
}

struct Complex {
    Interval re;
    Interval im;
};

Complex operator*(Complex a, Complex b) {
    return {a.re*b.re + -(a.im*b.im), a.re*b.im + a.im*b.re};
}

Complex exp_unit(Interval angle) {
    Complex out{{scale,scale}, {0,0}};
    Interval term{scale,scale};
    for (Int n=1; n<32; ++n) {
        term = term*angle;
        term = {floor_div(term.lo,n),ceil_div(term.hi,n)};
        switch (n % 4) {
            case 0: out.re = out.re + term; break;
            case 1: out.im = out.im + term; break;
            case 2: out.re = out.re + -term; break;
            case 3: out.im = out.im + -term; break;
        }
    }
    // prepare.py checks that the Taylor remainder is less than one unit.
    return {{out.re.lo-1,out.re.hi+1},{out.im.lo-1,out.im.hi+1}};
}

Int magnitude_up(Complex z) {
    const Int re = std::max(std::abs(z.re.lo), std::abs(z.re.hi));
    const Int im = std::max(std::abs(z.im.lo), std::abs(z.im.hi));
    const Wide square = Wide(re)*re + Wide(im)*im;
    if (square == 0) return 0;
    const auto high = static_cast<std::uint64_t>(square >> 64);
    const int bits = high ? 64+std::bit_width(high)
                          : std::bit_width(static_cast<std::uint64_t>(square));
    assert(bits <= 122);
    Int root = Int{1} << ((bits+1)/2);
    // Integer Newton iteration from above, followed by upward rounding.
    while (true) {
        const Int next = static_cast<Int>((root+square/root)/2);
        if (next >= root) return root+(Wide(root)*root != square);
        root = next;
    }
}

struct Prefix {
    std::size_t offset;
    std::string digits;
    std::array<std::vector<unsigned>,width> returns;
    std::array<Int,width> curvature{};
    // Each budget is followed by its nonnegative weights, in return order.
    std::vector<std::array<Int,width+1>> rows;
};

template <typename T>
void read(T& x) {
    if (!(std::cin >> x)) {
        std::cerr << "Incomplete certificate input\n";
        std::exit(2);
    }
}

void return_words(Prefix& prefix, const std::array<int,width>& critical,
                  std::size_t time, std::size_t state, unsigned mask) {
    if (state == 0) prefix.returns[time-1].push_back(mask);
    if (time == prefix.digits.size()) return;
    for (int digit=0; digit<=critical[state]; ++digit) {
        const std::size_t next = digit < critical[state] ? 0 : state+1;
        return_words(prefix, critical, time+1, next,
                     mask | (static_cast<unsigned>(digit) << (prefix.offset+time)));
    }
}

void initialise(Prefix& prefix, const std::array<int,width>& critical,
                const std::array<Int,width+1>& q, Int two_pi_up) {
    unsigned upper_path = 0;
    const std::size_t n = prefix.digits.size();
    for (std::size_t i=0; i<n; ++i) {
        if (prefix.digits[i] == '1') {
            return_words(prefix, critical, i+1, 0, upper_path);
            upper_path |= 1u << (prefix.offset+i);
        }
    }
    for (std::size_t r=0; r<n; ++r) {
        Int sum = 0;
        for (unsigned mask: prefix.returns[r]) {
            Int derivative = 0;
            for (std::size_t j=0; j<width; ++j) {
                if (mask & (1u<<j)) derivative += q[j];
            }
            const Int slope = mul_up(derivative,two_pi_up);
            sum = narrow(Wide(sum)+mul_up(slope,slope));
        }
        prefix.curvature[r] = mul_up(q[r+1],sum);
    }
}

std::array<Complex,masks> phases(Interval two_pi) {
    std::array<Complex,width> sites;
    for (auto& site: sites) {
        Int x;
        read(x);
        assert(-scale/2 <= x && x <= scale/2);
        site = exp_unit(Interval{x,x+1}*two_pi);
    }
    std::array<Complex,masks> all;
    all[0] = {{scale,scale},{0,0}};
    for (unsigned mask=1; mask<masks; ++mask) {
        const auto bit = static_cast<std::size_t>(std::countr_zero(mask));
        all[mask] = all[mask&(mask-1)]*sites[bit];
    }
    return all;
}

using Bounds = std::vector<std::array<Int,width>>;

Bounds norm_bounds(const std::vector<Prefix>& prefixes, std::array<Complex,masks> all,
                   const std::array<Int,width+1>& q) {
    Bounds out(prefixes.size());
    for (std::size_t i=0; i<prefixes.size(); ++i) {
        const auto& prefix = prefixes[i];
        for (std::size_t r=0; r<prefix.digits.size(); ++r) {
            Complex sum{};
            for (unsigned mask: prefix.returns[r]) {
                sum = {sum.re+all[mask].re,sum.im+all[mask].im};
            }
            out[i][r] = mul_up(q[r+1],magnitude_up(sum));
        }
    }
    return out;
}

int main(int argc, char** argv) {
    const bool pattern_margins = argc == 2 && std::string(argv[1]) == "--pattern-margins";
    if (argc != 1 && !pattern_margins) {
        std::cerr << "usage: check [--pattern-margins]\n";
        return 2;
    }
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    Interval two_pi;
    read(two_pi.lo);
    read(two_pi.hi);
    std::array<int,width> critical;
    for (int& d: critical) read(d);
    std::array<Int,width+1> q;
    for (Int& v: q) read(v);
    std::size_t total, gaps;
    read(total);
    read(gaps);
    std::vector<Prefix> prefixes(total);
    for (auto& prefix: prefixes) {
        read(prefix.offset);
        read(prefix.digits);
        assert(prefix.offset+prefix.digits.size() == width);
        initialise(prefix,critical,q,two_pi.hi);
        prefix.rows.resize(gaps);
        for (auto& row: prefix.rows) {
            read(row[0]);
            for (std::size_t r=0; r<prefix.digits.size(); ++r) read(row[r+1]);
        }
    }
    int patterns;
    read(patterns);
    int expected_patterns = 2;
    for (std::size_t i=2; i<width; ++i) expected_patterns *= 3;
    assert(patterns == expected_patterns);
    Int margin = std::numeric_limits<Int>::max();
    long long intervals = 0;
    for (int pat=0; pat<patterns; ++pat) {
        Int pattern_margin = std::numeric_limits<Int>::max();
        int pieces;
        read(pieces);
        for (int piece=0; piece<pieces; ++piece) {
            int count;
            Int mesh;
            read(count);
            read(mesh);
            const Int remainder_factor = ceil_div(mul_up(mesh,mesh),8);
            auto left = norm_bounds(prefixes,phases(two_pi),q);
            for (int interval=0; interval<count; ++interval) {
                auto right = norm_bounds(prefixes,phases(two_pi),q);
                for (std::size_t i=0; i<prefixes.size(); ++i) {
                    const auto& prefix = prefixes[i];
                    std::array<Int,width> remainder;
                    const std::size_t n = prefix.digits.size();
                    for (std::size_t r=0; r<n; ++r) {
                        remainder[r] = mul_up(remainder_factor,prefix.curvature[r]);
                    }
                    for (std::size_t g=0; g<gaps; ++g) {
                        const auto& row = prefix.rows[g];
                        Wide sum_left = 0, sum_right = 0, error = 0;
                        for (std::size_t r=0; r<n; ++r) {
                            const Int weight = row[n-r];
                            sum_left += Wide(left[i][r])*weight;
                            sum_right += Wide(right[i][r])*weight;
                            error += Wide(remainder[r])*weight;
                        }
                        const Int slack = row[0]-ceil_div(std::max(sum_left,sum_right)+error,scale);
                        pattern_margin = std::min(pattern_margin,slack);
                    }
                }
                left.swap(right);
                ++intervals;
            }
        }
        margin = std::min(margin,pattern_margin);
        if (pattern_margins) std::cout << "pattern " << pat << " " << pattern_margin << "\n";
    }
    std::string extra;
    assert(!(std::cin >> extra));
    std::cout << "patterns " << patterns << ", intervals " << intervals
              << ", rational window slack >= " << margin << "/" << scale << "\n";
    if (margin <= 0) return 1;
    std::cout << "CERTIFIED\n";
}
