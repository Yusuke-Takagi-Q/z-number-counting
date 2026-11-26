// Second verifier for the width-seven protocol emitted by prepare.py.
#include "../independent_interval.hpp"

#include <array>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using namespace crosscheck;
constexpr int width = 7;
constexpr int masks = 1 << width;
using Row = std::array<I, width>;
constexpr Box one{{scale, scale}, {}};

template<class T> T read() {
    T value{};
    if (!(std::cin >> value)) throw std::runtime_error("incomplete input");
    return value;
}
I input_unit() { return checked(Z(read<I>())*lift); }

struct Budget {
    Z cap;
    Row weights{};
};
struct Prefix {
    int offset;
    std::array<std::vector<int>, width> endings;
    Row curvature{};
    std::vector<Budget> budgets;
};

bool returns(const std::string& top, const std::array<int,width>& critical,
             int word, int length) {
    int age = -1; // still on the upper path
    for (int j = 0; j < length; ++j) {
        const int digit = (word >> j) & 1;
        const int limit = age < 0 ? top[j]-'0' : critical[age];
        if (digit > limit) return false;
        if (digit < limit) age = 0;
        else if (age >= 0) ++age;
    }
    return age == 0;
}

int main() try {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    const Range pi2{input_unit(), input_unit()};
    if (pi2.lo <= 0 || pi2.lo > pi2.hi || pi2.hi >= 8*scale)
        throw std::runtime_error("invalid pi enclosure");
    std::array<int,width> critical;
    for (auto& d : critical) {
        d = read<int>();
        if (d != 0 && d != 1) throw std::runtime_error("invalid critical digit");
    }
    std::array<I,width+1> q;
    for (auto& v : q) {
        v = input_unit();
        if (v <= 0 || v > scale) throw std::runtime_error("invalid geometric weight");
    }

    const int count = read<int>(), gaps = read<int>();
    if (count <= 0 || gaps <= 0) throw std::runtime_error("invalid budget dimensions");
    std::vector<Prefix> prefixes;
    for (int p = 0; p < count; ++p) {
        Prefix prefix{};
        prefix.offset = read<int>();
        const auto top = read<std::string>();
        const int length = static_cast<int>(top.size());
        if (prefix.offset < 0 || length < 1 || length > width ||
            prefix.offset != width-length || top.find_first_not_of("01") != std::string::npos)
            throw std::runtime_error("invalid prefix");
        for (int r = 1; r <= length; ++r) {
            // Enumerate digit words directly and simulate each one.
            Z derivative_squares = 0;
            for (int word = 0; word < (1<<r); ++word) {
                if (!returns(top,critical,word,r)) continue;
                const int mask = word << prefix.offset;
                prefix.endings[r-1].push_back(mask);
                I derivative = 0;
                for (int j = 0; j < width; ++j)
                    if ((mask >> j) & 1) derivative = checked(Z(derivative)+q[j]);
                const I slope = up(Z(derivative)*pi2.hi,scale);
                derivative_squares += up(Z(slope)*slope,scale);
            }
            prefix.curvature[r-1] = up(derivative_squares*q[r],scale);
        }
        for (int g = 0; g < gaps; ++g) {
            Budget budget{Z(input_unit())*scale, {}};
            for (int r = length-1; r >= 0; --r) {
                budget.weights[r] = input_unit();
                if (budget.weights[r] < 0 || budget.weights[r] > 8*scale)
                    throw std::runtime_error("invalid renewal weight");
            }
            prefix.budgets.push_back(budget);
        }
        prefixes.push_back(std::move(prefix));
    }

    auto read_endpoint = [&]() {
        std::array<Box,width> twists;
        for (auto& twist : twists) {
            const I lo = input_unit();
            if (lo < -scale/2 || lo > scale/2)
                throw std::runtime_error("phase outside centered interval");
            twist = horner_exponential(multiply({lo,checked(Z(lo)+lift)},pi2));
        }
        // Build subset products by successively adjoining the next digit.
        std::array<Box,masks> products;
        products[0] = one;
        for (int j = 0; j < width; ++j)
            for (int m = 0; m < (1<<j); ++m)
                products[m+(1<<j)] = products[m]*twists[j];

        std::vector<Row> result(prefixes.size());
        for (std::size_t p = 0; p < prefixes.size(); ++p) {
            for (int r = 0; r < width-prefixes[p].offset; ++r) {
                Box sum{};
                for (const int mask : prefixes[p].endings[r]) sum = sum + products[mask];
                result[p][r] = up(Z(modulus_up(sum))*q[r+1],scale);
            }
        }
        return result;
    };

    I minimum = INT64_MAX;
    I intervals = 0;
    const int patterns = read<int>();
    if (patterns != 486) throw std::runtime_error("incorrect number of patterns");
    for (int pattern = 0; pattern < patterns; ++pattern) {
        const int pieces = read<int>();
        if (pieces <= 0 || pieces > width) throw std::runtime_error("invalid phase partition");
        for (int piece = 0; piece < pieces; ++piece) {
            const int steps = read<int>();
            const I mesh = input_unit();
            if (steps <= 0 || mesh <= 0 || mesh > scale)
                throw std::runtime_error("invalid mesh");
            const I chord = up(Z(mesh)*mesh,8*Z(scale));
            std::vector<std::vector<Z>> errors(prefixes.size());
            for (std::size_t p = 0; p < prefixes.size(); ++p) {
                for (const auto& budget : prefixes[p].budgets) {
                    Z correction = 0;
                    for (int r = 0; r < width-prefixes[p].offset; ++r)
                        correction += Z(budget.weights[r])*up(Z(chord)*prefixes[p].curvature[r],scale);
                    errors[p].push_back(correction);
                }
            }
            auto left = read_endpoint();
            for (int j = 0; j < steps; ++j) {
                auto right = read_endpoint();
                for (std::size_t p = 0; p < prefixes.size(); ++p) {
                    for (int g = 0; g < gaps; ++g) {
                        const auto& budget = prefixes[p].budgets[g];
                        Z lhs_left = 0, lhs_right = 0;
                        for (int r = 0; r < width-prefixes[p].offset; ++r) {
                            lhs_left += Z(budget.weights[r])*left[p][r];
                            lhs_right += Z(budget.weights[r])*right[p][r];
                        }
                        const Z slack = budget.cap-std::max(lhs_left,lhs_right)-errors[p][g];
                        minimum = std::min(minimum,down(slack,scale));
                    }
                }
                left.swap(right);
                ++intervals;
            }
        }
    }
    std::string extra;
    if (std::cin >> extra) throw std::runtime_error("unexpected trailing input");
    std::cout << patterns << " patterns, " << intervals << " intervals; minimum slack >= "
              << minimum << "/" << scale << "\n";
    if (minimum <= 0) return 1;
    std::cout << "INDEPENDENT WINDOW CHECK PASSED\n";
} catch (const std::exception& e) {
    std::cerr << e.what() << "\n";
    return 2;
}
