"""Compute a rational interval bound from a verified tile certificate."""
import argparse
import json
from fractions import Fraction as F
from pathlib import Path


def interval_bound(cert, delta, *, terminal_block=False):
    q = F(2, 3)
    B = F(cert.get("B", "631/1000"))
    T, M, D = (cert[key] for key in ("T", "M", "D"))
    # The original T10 certificate predates cell-count metadata.
    if "cells" in cert:
        cells = cert["cells"]
    elif (T, M, D, cert["mesh"]) == (10, 9, 18, 7):
        cells = 128
    else:
        raise ValueError("the verified number of phase cells is required")
    Lambda, rho, kappa, log_upper = (
        F(cert[key]) for key in ("Lambda", "rho", "kappa", "S_upper")
    )
    c = [[F(value) for value in row] for row in cert["c"]]
    cmax = max(F(1), max(map(max, c)))
    a_upper = F(117, 200)
    assert 3**200 < 2**317
    assert B >= F(1, 2) and B < 1
    exponent = 1 / B
    assert 2**exponent.numerator < 3**exponent.denominator
    assert q < kappa < 1 and rho >= 1 and q**T * rho < 1
    assert Lambda**log_upper.denominator <= 2**(T * log_upper.numerator)
    assert (F(3, 2))**T <= Lambda < (F(11, 7))**T
    # These bounds imply lambda^3 < 4 and lambda^B > 5/4.
    assert F(11, 7)**3 < 4 and B >= F(3, 5)
    assert F(3, 2)**3 > F(5, 4)**5

    # Sufficient tail tests, using the number of frequencies 2^d.
    for exponent in (M * (1 - B), M * B):
        assert 2**exponent.denominator <= F(3, 2)**exponent.numerator

    # Compute c(y_a) for a finite initial segment. The geometric tail
    # uses y_a c(y_a) <= cmax, which follows from q^T rho < 1.
    limit = max(64, T)
    y = [F(1)]
    digits = []
    for _ in range(limit + D + 128):
        value = y[-1] / q
        digit = value // 1
        digits.append(digit)
        y.append(value - digit)
    lookup = {word: i for i, word in enumerate(cert["cylinders"])}

    def potential_at_age(age):
        if age == 0:
            return F(1)
        first_one = digits.index(1, age)
        word = "".join(map(str, digits[first_one:first_one + D]))
        quotient, remainder = divmod(first_one - age, T)
        return rho**quotient * c[remainder][lookup[word]]

    ratio = q / kappa
    weights = [ratio**age * y[age] for age in range(limit)]
    tail = ratio**limit / (1 - ratio)
    G = sum(weights) + tail
    Z = sum(weight * potential_at_age(age) for age, weight in enumerate(weights)) + cmax * tail
    tail_ratio = 3**T * kappa**(M * T) / Lambda
    assert tail_ratio < 1
    C1 = T * 3**(T - 1) / (1 - tail_ratio)**2 * (
        1 + cells * sum((Z / (weights[j] * potential_at_age(j)))**M for j in range(T))
    )
    # The r truncated blocks contribute an extra factor C1**(-r/M).
    # Bound its sum by a geometric series when C1 >= 100**M.
    truncated = F(1, 99) if C1 >= 100**M else M
    block = 3 * (1 + truncated * q**-2) * G**(M + 1) * C1
    small_depth = -(-2 // (1 - B))
    fourier = 20 * block + 2**small_depth - 2
    searches = M - 1 if terminal_block else M
    theta = (1 - B)**searches / (2 - B)
    eta = theta / (1 - theta)
    saving = (1 - a_upper - B * log_upper) * eta
    assert 0 < delta <= saving
    count = 45 * (1 + 2 * fourier) + 1
    return {"G": G, "Z": Z, "C1": C1, "fourier": fourier, "eta": eta,
            "saving_lower": saving, "count_upper": count}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificate", type=Path)
    parser.add_argument("--delta", default="0.0000023922", type=F)
    parser.add_argument("--constant", type=int)
    parser.add_argument("--terminal-block", action="store_true",
                        help="use the terminal injective block from note 19")
    args = parser.parse_args()
    cert = json.loads(args.certificate.read_text())
    bounds = interval_bound(cert, args.delta, terminal_block=args.terminal_block)
    if args.constant is not None:
        assert bounds["count_upper"] <= args.constant
    for key, value in bounds.items():
        print(f"{key} ~ {float(value):.8g}")
