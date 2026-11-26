"""Generate rationally enclosed input for check.cpp, for the 3/2 age potential."""
import argparse
from fractions import Fraction as F
from itertools import product
import json
from math import ceil, factorial, floor
from pathlib import Path
import sys

W, DEPTH, CYLINDER_DEPTH, G, N_PLUS = 7, 8, 18, 48, 30
INV_BETA = F(2, 3)
SCALE = 2**42
PARAMS = Path(__file__).resolve().parents[1] / "params.json"


def lower(value):
    return floor(SCALE * value)


def upper(value):
    return ceil(SCALE * value)


def atan_bounds(x, terms):
    partial = sum((-1)**j * x**(2*j+1) / (2*j+1) for j in range(terms))
    following = partial + (-1)**terms * x**(2*terms+1)/(2*terms+1)
    return min(partial, following), max(partial, following)


def pi_bounds():
    alo, ahi = atan_bounds(F(1, 5), 18)
    blo, bhi = atan_bounds(F(1, 239), 8)
    return 16*alo-4*bhi, 16*ahi-4*blo


class Certificate:
    def __init__(self, params, rho, tail_rho):
        self.rho = F(params["rho"] if rho is None else rho)
        self.tail_rho = self.rho if tail_rho is None else F(tail_rho)
        rho = self.rho
        self.cap = F(params["K0"])
        self.kappa = {tuple(map(int, key)): F(value)
                      for key, value in params["kappa"].items()}
        self.lam = {int(key): F(value) for key, value in params["lambda"].items()}
        self.lam[0] = max(self.kappa.values())
        assert INV_BETA**W < self.tail_rho <= rho < 1
        assert all(0 < x <= self.cap for x in self.kappa.values())
        assert all(0 < self.lam[i] <= self.lam[i+1] for i in range(1, W-1))
        self.fill_weights(W)
        assert self.lam[W-1] <= self.lam[W]
        assert min(F(1), min(self.kappa.values()), self.lam[1]) > F(94, 100)
        assert self.cap < F(139, 100) and rho * F(103, 100)**W > 1

        self.digits = []
        self.y = [F(1)]
        for j in range(1, 900):
            value = 3*self.y[-1]/2
            self.digits.append(floor(value))
            self.y.append(value-self.digits[-1])
        self.q = [INV_BETA**a for a in range(900)]
        self.f = [F(0)] + [d*self.q[j+1] for j, d in enumerate(self.digits)]
        self.tail = [power*y for power, y in zip(self.q, self.y)]
        self.cylinders = self.words(CYLINDER_DEPTH)
        self.bound_wc = self.sup_wc((), 0, False)

        # u(n) is an integer return count multiplied by q^n.
        self.prepare_renewals()
        self.pc = [[self.end_weight(a, g) for a in range(W+1)]
                   for g in range(G+2)]

    def fill_weights(self, m):
        for t in range(len(self.lam), m+1):
            self.lam[t] = max(self.lam[t-W], self.cap)/self.tail_rho

    def leading_zeros(self, m):
        self.fill_weights(m)
        return self.lam[m]

    def words(self, length):
        level = [((), 0, F(0))]
        for i in range(length):
            level = [(word+(digit,), 0 if digit < self.digits[state] else state+1,
                      left+digit*self.q[i+1])
                     for word, state, left in level
                     for digit in range(self.digits[state]+1)]
        return [row for row in level if row[0][0] == 1]

    def state(self, word, start=0):
        for digit in word:
            assert digit <= self.digits[start]
            start = 0 if digit < self.digits[start] else start+1
        return start

    def sup_wc(self, suffix, start, with_cap):
        """Upper bound for w*c(w), or w*max(c(w), cap), on the cylinder."""
        length = len(suffix)
        if 1 in suffix:
            first = suffix.index(1)
            if first:
                weight = self.leading_zeros(first)
            else:
                key = suffix[:DEPTH]
                weight = max(val for word, val in self.kappa.items()
                             if word[:len(key)] == key)
            if with_cap:
                weight = max(weight, self.cap)
            left = sum(d*self.q[i+1] for i, d in enumerate(suffix))
            right = left+self.q[length]*self.y[self.state(suffix, start)]
            return right*weight

        # For m>=W, c_{m+W}=c_m/tail_rho, so q^m c_m decreases along each residue class.
        last = max(length, W)+W
        return max(self.q[m]*max(self.leading_zeros(m), self.cap if with_cap else 0)
                   for m in range(length, last))

    def orbit_weight(self, a):
        if a == 0:
            return F(1)
        first = self.digits.index(1, a)-a
        if first:
            return self.leading_zeros(first)
        return self.kappa[tuple(self.digits[a:a+DEPTH])]

    def prepare_renewals(self):
        reference, band, cut = 560, 128, 130
        end = reference+cut
        counts = [1]
        for n in range(1, end):
            counts.append(sum(counts[n-j] for j in range(1, n+1) if self.digits[j-1]))
        u = [count*power for count, power in zip(counts, self.q)]
        self.orbit = [self.orbit_weight(a) for a in range(end)]
        coefficients = [T*c for T, c in zip(self.tail, self.orbit)]
        self.psi = [sum(u[n-a]*coefficients[a] for a in range(n+1)) for n in range(end)]
        u_far = max(u[reference-band+1:reference+1])+self.q[band]/(1-INV_BETA)
        psi_far = (u_far*sum(coefficients[:cut+1])
                   + self.q[cut+1]*self.bound_wc/(1-INV_BETA))
        self.psi_max = max(max(self.psi), psi_far)
        self.psi_plus = max(max(self.psi[N_PLUS:]), psi_far)
        assert self.psi_max <= self.cap

    def escapes(self, word, gap):
        last = min(len(word), gap) if gap <= G else len(word)
        if gap <= G:
            return sum(self.q[i+1]*self.psi[gap-1-i] for i in range(last) if word[i])
        assert G >= len(word)+N_PLUS-1
        return self.psi_plus*sum(self.q[i+1] for i in range(last) if word[i])

    def top(self, word, prefix, gap):
        """Unnormalised upper-path term, with prefix already read."""
        suffix = word[prefix:]
        if gap <= min(len(suffix), G):
            residual = suffix[gap:]
            state = self.state(word[:prefix+gap])
            rest = self.q[gap]*self.sup_wc(residual, state, False)
        else:
            state = self.state(word)
            rest = self.q[len(suffix)]*self.sup_wc((), state, True)
        return self.q[prefix]*(self.escapes(suffix, gap)+rest)

    def end_weight(self, a, g):
        """Upper bound on T_a (P^g c)(y_a), with no normalisation."""
        if g <= G:
            return (sum(self.f[a+i]*self.psi[g-i] for i in range(1, g+1))
                    + self.tail[a+g]*self.orbit[a+g])
        split = G-N_PLUS
        early = sum(self.f[a+i] for i in range(1, split+1))
        return (early*self.psi_plus+(self.tail[a]-early)*self.psi_max
                + self.q[a+G+1]*self.bound_wc)

    def budgets(self):
        groups = {}
        cap_check = max(self.psi_max, max(self.kappa.values()))
        for word, _, left in self.cylinders:
            for g in range(1, G+2):
                cap_check = max(cap_check, self.top(word, 0, g)/left)
            for m in range(W):
                size = W-m
                prefix = word[:size]
                weight = self.kappa[word[:DEPTH]] if m == 0 else self.lam[m]
                key = (m, prefix)
                current = [self.rho*weight*left-self.top(word, size, g)
                           for g in range(G+2)]
                if key in groups:
                    current = list(map(min, zip(current, groups[key])))
                groups[key] = current
        assert cap_check <= self.cap
        print("rational cap slack >=", lower(self.cap-cap_check), "/", SCALE, file=sys.stderr)

        # The initial state s=1 uses the same G as its critical prefix.
        key = (0, tuple(self.digits[:W]))
        critical = [self.rho-self.pc[g][W] for g in range(G+2)]
        groups[key] = list(map(min, zip(critical, groups.get(key, critical))))
        return groups


def pieces(pattern):
    """Affine phase vectors on the closed one-sided pieces of [-1/2,1/2)."""
    parts = [(F(-1, 2), F(1, 2), [(F(1), F(0))])]
    for r in pattern:
        new = []
        for a, b, coefficients in parts:
            slope, shift = coefficients[-1]
            slope, shift = INV_BETA*slope, INV_BETA*shift-F(r, 3)
            cuts = sorted({a, b} | {
                (F(n)+F(1, 2)-shift)/slope
                for n in (-2, -1, 0)
                if a < (F(n)+F(1, 2)-shift)/slope < b
            })
            for left, right in zip(cuts, cuts[1:]):
                offset = floor(slope*(left+right)/2+shift+F(1, 2))
                new.append((left, right, coefficients+[(slope, shift-offset)]))
        parts = new
    return parts


def output(args, params):
    cert = Certificate(params, args.rho, args.tail_rho)
    budgets = cert.budgets()
    pi_lo, pi_hi = pi_bounds()
    assert pi_hi < F(22, 7) and F(22, 7)**32/factorial(32) < F(1, SCALE)
    print(lower(2*pi_lo), upper(2*pi_hi))
    print(*cert.digits[:W])
    print(*(upper(INV_BETA**r) for r in range(W+1)))
    print(len(budgets), G+2)
    for (m, prefix), rhs in sorted(budgets.items()):
        print(m, "".join(map(str,prefix)))
        for g, bound in enumerate(rhs):
            print(lower(bound), *(upper(cert.pc[g][a]) for a in range(len(prefix))))
    patterns = product(*(range(1, 3) if j == args.pivot else range(3)
                         for j in range(W-1)))
    print(2*3**(W-2))
    for pattern in patterns:
        parts = pieces(pattern)
        print(len(parts))
        for left, right, affine in parts:
            count = ceil((right-left)*args.grid)
            width = (right-left)/count
            print(count, upper(width))
            for i in range(count+1):
                x = left+i*width
                phases = [-(slope*x+shift) for slope, shift in affine]
                assert all(-F(1,2) <= phase <= F(1,2) for phase in phases)
                print(*(lower(phase) for phase in phases))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--params", default=PARAMS, type=Path)
    parser.add_argument("--rho", help="override the value in params.json")
    parser.add_argument("--tail-rho", help="override the multiplier in the deep-age recurrence")
    parser.add_argument("--grid", default=480, type=int)
    parser.add_argument("--pivot", default=0, type=int)
    args = parser.parse_args()
    params = json.loads(args.params.read_text(), parse_float=F)
    W, DEPTH, CYLINDER_DEPTH, G, N_PLUS = (params[key] for key in ("W", "D", "DC", "G", "NPLUS"))
    if not 3 <= W <= 7 or not 0 <= args.pivot < W-1:
        parser.error("width must be in [3, 7] and pivot in [0, width-2]")
    if args.grid < 1:
        parser.error("--grid must be positive")
    output(args, params)
