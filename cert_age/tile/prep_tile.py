"""Rigorous check of the tile-moment condition (C) for base 3/2 (notes 15-16).

Rational preprocessing, a call to tilecheck (integer interval arithmetic), and the final
rational comparisons.  Standard library only.

usage: python3 prep_tile.py c_table.json --lam-exp 0.43 --mesh 6 --threads 8 --bin ./tilecheck
  c_table.json : output of cylopt (rounded to rationals) or a saved certificate (exact fractions)
  --lam-exp    : choose a rational Lambda near e^(T*lam_exp); all certification uses that rational
  --mesh k     : cells partition [-1/2,1/2): grid of width 2^-k plus cut points of depth <= --cut-depth
"""
import argparse, json, subprocess, sys, time
from fractions import Fraction as F
from math import comb, exp, floor, ceil, factorial

q = F(2, 3); SCALE = 2**42

def atan_bounds(x, terms):
    s = sum((-1)**j * x**(2*j+1) / (2*j+1) for j in range(terms))
    t = s + (-1)**terms * x**(2*terms+1) / (2*terms+1)
    return min(s, t), max(s, t)

def pi_bounds():
    alo, ahi = atan_bounds(F(1, 5), 18); blo, bhi = atan_bounds(F(1, 239), 8)
    return 16*alo - 4*bhi, 16*ahi - 4*blo

def refine_potential(ctab, P, E, NJ, plist, classes, cyls, Lg, Rg, dig, y, T, M, D, one_pid,
                     class_of_digits, iters, frozen):
    """Floating-point Perron iteration for c with the power sums P fixed (choice of c only)."""
    import bisect
    qf = 2 / 3; G = len(cyls)
    c = [[float(v) for v in row] for row in ctab]
    Pmax = {pp: [max(P[(j, pp)][k] for j in range(NJ)) / SCALE for k in range(M + 1)] for pp in range(len(plist))}
    Emax = max(E.values()); Lg = [float(v) for v in Lg]; Rg = [float(v) for v in Rg]
    binom = [comb(M, k) for k in range(M + 1)]
    def poly(Pv, N): return sum(binom[k] * N**(M - k) * Pv[k] for k in range(M + 1))
    m1, g1 = class_of_digits(dig[T:T+D+200]); yT = float(y[T])
    Lam = None
    for it in range(iters):
        rho = ((Emax / Lam) ** (1.0 / M)) if Lam else 1.6
        def cf(m, g):
            f = 1.0
            while m >= T: m -= T; f *= rho
            return f * c[m][g]
        maxc_all = [max(row) for row in c]
        def mx_m(m):
            f = 1.0
            while m >= T: m -= T; f *= rho
            return f * maxc_all[m]
        def mx_pre(m, pre):
            lo = bisect.bisect_left(cyls, pre); hi = bisect.bisect_left(cyls, pre + "9")
            f = 1.0
            while m >= T: m -= T; f *= rho
            return f * max(c[m][g] for g in range(lo, hi))
        need = {}; ub = 0
        for (m, g, pp) in classes:
            w = cyls[g]; al = w[T-m:]
            Lp = sum(int(w[i]) * qf**(i+1) for i in range(T - m)); zhi = (Rg[g] - Lp) / qf**(T - m)
            if "1" in al:
                mz = al.index("1"); N = qf**T * min(zhi, qf**mz) * mx_pre(mz, al[mz:])
            else:
                la = len(al); N = qf**T * max(qf**mm * mx_m(mm) for mm in range(la, la + T))
            v = poly(Pmax[pp], N); slo = qf**m * Lg[g]
            need[(m, g)] = v ** (1.0 / M) / slo; ub = max(ub, v / (slo * c[m][g])**M)
        v1 = poly(Pmax[one_pid], qf**T * yT * cf(m1, g1)); ub = max(ub, v1)
        Lam = ub; n1 = v1 ** (1.0 / M)
        for (m, g) in need:
            if (m, g) not in frozen:
                c[m][g] = (c[m][g] * need[(m, g)] / n1) ** 0.5
    import math
    print(f"  refine: ln(Lambda)/T <= {math.log(ub)/T:.4f} (with P_j fixed, Pmax over cells)", flush=True)
    return [[ctab[m][g] if (m, g) in frozen else F(round(v * 10**9), 10**9)
             for g, v in enumerate(row)] for m, row in enumerate(c)]

def main():
    if not __debug__:
        raise SystemExit("Run prep_tile.py with Python assertions enabled.")
    ap = argparse.ArgumentParser()
    ap.add_argument("table"); ap.add_argument("--lam-exp", type=float)
    ap.add_argument("--mesh", type=int); ap.add_argument("--cut-depth", type=int); ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--bin", default="./tilecheck"); ap.add_argument("--kappa")
    ap.add_argument("--no-symmetry", action="store_true",
                    help="compute every phase cell independently")
    ap.add_argument("--B"); ap.add_argument("--save", default=None)
    ap.add_argument("--refine", type=int, default=0); ap.add_argument("--cache", default=None);
    ap.add_argument("--freeze-returns", action="store_true",
                    help="keep c(y_a), a < T, fixed when refining to reuse the cached sums")
    ap.add_argument("--table-out", default=None); ap.add_argument("--refine-iters", type=int, default=40)
    args = ap.parse_args()
    tab = json.load(open(args.table)); T, M, D = tab["T"], tab["M"], tab["D"]
    replay = args.lam_exp is None
    if replay and "Lambda" not in tab:
        ap.error("--lam-exp is required for a potential table without a saved Lambda")
    B = F(args.B or tab.get("B", "631/1000"))
    kappa = F(args.kappa or tab.get("kappa", "4/5"))
    # potential as exact rationals (any positive rational table is a valid choice)
    cyls = tab["cylinders"]
    ctab = [[F(v) if isinstance(v, str) else F(round(v * 10**9), 10**9)
             for v in row] for row in tab["c"]]
    assert len(ctab) == T and all(len(row) == len(cyls) for row in ctab)
    assert all(v > 0 for row in ctab for v in row)
    Lam = F(tab["Lambda"]) if replay else F(ceil(exp(T * args.lam_exp) * 10**9), 10**9)
    assert Lam > 1 and q < kappa < 1

    # critical digits and y_a
    y = [F(1)]; dig = []
    for _ in range(T + D + 400):
        v = F(3, 2) * y[-1]; dig.append(floor(v)); y.append(v - dig[-1])
    def state_after(word, st=0):
        for x in word:
            assert x <= dig[st]; st = 0 if x < dig[st] else st + 1
        return st
    # cylinders: regenerate and compare with the table order
    level = [((), 0)]
    for _ in range(D):
        level = [(w + (x,), 0 if x < dig[st] else st + 1) for w, st in level for x in range(dig[st] + 1)]
    words = sorted("".join(map(str, w)) for w, st in level if w[0] == 1)
    assert words == cyls, "cylinder list mismatch"
    G = len(cyls); index = {w: g for g, w in enumerate(cyls)}
    Lg = [sum(int(w[i]) * q**(i+1) for i in range(D)) for w in cyls]
    Rg = [Lg[g] + q**D * y[state_after(tuple(map(int, cyls[g])))] for g in range(G)]

    def class_of_digits(ds):          # ds: greedy digits of a number in (0,1)
        m = 0
        while ds[m] == 0: m += 1
        return m, index["".join(map(str, ds[m:m+D]))]

    # c(y_a) for a<T (y_0 = 1 has c = 1)
    cy = [F(1)]
    for a in range(1, T + 1):
        m, g = class_of_digits(dig[a:a+D+200]); assert m < T, "c(y_a) needs m<T"
        cy.append(ctab[m][g])
    wgt = [ceil(SCALE * q**T * y[a] * cy[a]) for a in range(T)]

    pi_lo, pi_hi = pi_bounds()
    assert pi_hi < F(22, 7) and F(22, 7)**32 / factorial(32) < F(1, SCALE)
    tw_lo, tw_hi = floor(SCALE * 2 * pi_lo), ceil(SCALE * 2 * pi_hi)

    # prefixes: classes (m,g), m<T, plus s=1
    prefixes = {}
    def pid(p):
        if p not in prefixes: prefixes[p] = len(prefixes)
        return prefixes[p]
    classes = []
    for m in range(T):
        for g in range(G):
            p = "0" * m + cyls[g][:T-m]
            classes.append((m, g, pid(p)))
    one_pid = pid("".join(map(str, dig[:T])))
    plist = sorted(prefixes, key=prefixes.get)

    # cells: uniform grid of width 2^-k together with all cut points of the phase maps of depth <= i0
    k = tab.get("mesh", 6) if args.mesh is None else args.mesh
    i0 = tab.get("cut_depth", 4) if args.cut_depth is None else args.cut_depth
    KX = max(k, T + 2, i0 + 2)
    pts = {F(-1, 2) + F(i, 2**k) for i in range(2**k + 1)}
    pieces = [(F(-1, 2), F(1, 2), F(1), F(0))]           # x-interval, phi = A x + Bc
    for depth in range(1, i0 + 1):
        nxt = []
        for a, b, A_, B_ in pieces:
            for r in (0, 1, 2):
                lo_v, hi_v = A_ * a + B_, A_ * b + B_
                cuts = []
                if r == 1 and lo_v < F(1, 4) < hi_v: cuts = [(F(1, 4) - B_) / A_]
                if r == 2 and lo_v < F(-1, 4) < hi_v: cuts = [(F(-1, 4) - B_) / A_]
                ends = [a] + cuts + [b]
                for u, v in zip(ends, ends[1:]):
                    mid = A_ * (u + v) / 2 + B_
                    cc = 0 if r == 0 else ((1 if mid < F(1, 4) else -2) if r == 1 else (-1 if mid >= F(-1, 4) else 2))
                    nxt.append((u, v, 2 * A_ / 3, (2 * B_ + cc) / 3)); pts.update([u, v])
        pieces = nxt
    pts = sorted(pts)
    assert all((x * 2**KX).denominator == 1 for x in pts) and pts[0] == F(-1, 2) and pts[-1] == F(1, 2)
    bnd = [int(x * 2**KX) for x in pts]; Js = list(zip(bnd, bnd[1:]))
    assert F(2, 3)**T <= F(1, 2)

    def curvature(pstr):
        # Record curvature at each return time. A path ending at age a
        # has the same modulus as its return prefix at time T-a.
        Z = (0, F(0), F(0)); ages = [Z] * T; top = (1, F(0), F(0))
        returns = []
        def add(u, v): return (u[0]+v[0], u[1]+v[1], u[2]+v[2])
        def one(u, i):   # append a digit 1 at step i
            n, s1, s2 = u; e = q**i
            return (n, s1 + n*e, s2 + 2*e*s1 + n*e*e)
        for i in range(1, T + 1):
            new = [Z] * T
            for a in range(T - 1):
                if ages[a][0] == 0: continue
                if dig[a] == 1: new[0] = add(new[0], ages[a]); new[a+1] = add(new[a+1], one(ages[a], i))
                else: new[a+1] = add(new[a+1], ages[a])
            if pstr[i-1] == "1": new[0] = add(new[0], top); top = one(top, i)
            returns.append(new[0])
            ages = new
        return [ceil(SCALE * (2*pi_hi)**2 * u[2]) for u in reversed(returns)]
    curv_lines = [" ".join(map(str, curvature(s))) for s in plist]
    inp = [f"{T} {M} {KX} {args.threads}", f"{tw_lo} {tw_hi}", " ".join(map(str, wgt)),
           " ".join(map(str, dig[:T])), str(len(bnd)), " ".join(map(str, bnd))] + \
          [str(len(plist))] + plist + curv_lines
    print(f"T={T} M={M} D={D} cylinders={G} classes={len(classes)} prefixes={len(plist)} "
          f"mesh=2^-{k}+cuts(depth {i0}) ({len(Js)} cells)  Lambda={Lam}", flush=True)

    def validate_sums(E, P):
        assert set(E) == set(range(len(Js))) and all(e >= 3**T for e in E.values())
        assert len(P) == len(Js) * len(plist)
        for j, e in E.items():
            for p in range(len(plist)):
                row = P[j, p]
                assert len(row) == M + 1 and row[0] == e * SCALE
                assert all(v >= 0 for v in row)
        return E, P

    def run_tilecheck(wgt):
        inp[2] = " ".join(map(str, wgt))
        checker = [args.bin]
        if args.no_symmetry:
            checker.append("--no-symmetry")
        key = "\n".join(inp) + ("\n--no-symmetry" if args.no_symmetry else "")
        if args.cache:
            try:
                ch = json.load(open(args.cache))
                if ch["key"] == key:
                    print("using cached tilecheck output", flush=True)
                    E = {int(k): v for k, v in ch["E"].items()}
                    P = {tuple(map(int, k.split(","))): v for k, v in ch["P"].items()}
                    return validate_sums(E, P)
            except FileNotFoundError: pass
        t0 = time.time()
        out = subprocess.run(checker, input="\n".join(inp) + "\n", capture_output=True, text=True)
        assert out.returncode == 0, out.stderr[-2000:]
        print(f"tilecheck finished in {time.time()-t0:.0f}s", flush=True)
        E = {}; P = {}
        for line in out.stdout.split("\n"):
            if not line: continue
            f = line.split()
            if f[0] == "E": E[int(f[1])] = int(f[2])
            else: P[(int(f[1]), int(f[2]))] = [int(v) for v in f[3:]]
        validate_sums(E, P)
        if args.cache:
            json.dump({"key": key, "E": E, "P": {f"{a},{b}": v for (a, b), v in P.items()}}, open(args.cache, "w"))
        return E, P
    E, P = run_tilecheck(wgt)
    frozen = ({class_of_digits(dig[a:a+D+200]) for a in range(1, T)}
              if args.freeze_returns else set())
    for outer in range(args.refine):
        # Perron updates of the potential (floating point) with the rigorous power sums fixed.
        ctab = refine_potential(ctab, P, E, len(Js), plist, classes, cyls, Lg, Rg, dig, y, T, M, D,
                                one_pid, class_of_digits, args.refine_iters, frozen)
        cy = [F(1)] + [ctab[m][g] for (m, g) in (class_of_digits(dig[a:a+D+200]) for a in range(1, T + 1))]
        wgt = [ceil(SCALE * q**T * y[a] * cy[a]) for a in range(T)]
        E, P = run_tilecheck(wgt)
    if args.table_out:
        json.dump({"T": T, "M": M, "D": D, "cylinders": cyls, "c": [[float(v) for v in row] for row in ctab]},
                  open(args.table_out, "w"))
    Emax = max(E.values())
    # Choose rho on the 10^-9 grid and check condition (R) exactly.
    step = F(1, 10**9)
    if replay:
        rho = F(tab["rho"])
    else:
        rho = F(ceil((float(Emax) / float(Lam))**(1.0 / M) * 10**9), 10**9)
        while rho**M * Lam < Emax:
            rho += step
        while rho > step and (rho - step)**M * Lam >= Emax:
            rho -= step
    slack = rho**M * Lam - Emax
    assert rho > 0 and slack >= 0 and q**T * rho < 1
    print(f"pieces: max_J E(J) = {Emax}; rho = {rho}; condition (R) slack ~ {float(slack):.9g}")
    if rho > step and (rho / step).denominator == 1:
        previous_slack = (rho - step)**M * Lam - Emax
        print(f"Condition (R) slack at rho - 10^-9: {float(previous_slack):.9g}")

    def cfun(m, g):
        f = F(1)
        while m >= T: m -= T; f *= rho
        return f * ctab[m][g]
    maxc_all = [max(ctab[m]) for m in range(T)]
    def maxc_m(m):
        f = F(1)
        while m >= T: m -= T; f *= rho
        return f * maxc_all[m]
    # ranges of cylinders with a given prefix (cylinders are sorted)
    import bisect
    def maxc_prefix(m, pre):
        lo = bisect.bisect_left(cyls, pre); hi = bisect.bisect_left(cyls, pre + "9")
        assert hi > lo
        f = F(1); mm = m
        while mm >= T: mm -= T; f *= rho
        return f * max(ctab[mm][g] for g in range(lo, hi))

    def N_bound(m, g):
        w = cyls[g]; alpha = w[T-m:]
        Lp = sum(int(w[i]) * q**(i+1) for i in range(T - m))
        zhi = (Rg[g] - Lp) / q**(T - m)
        if "1" in alpha:
            mz = alpha.index("1")
            return q**T * min(zhi, q**mz) * maxc_prefix(mz, alpha[mz:])
        la = len(alpha)
        return q**T * max(q**mm * maxc_m(mm) for mm in range(la, la + T))

    # judge
    Pmax = {}
    for p in range(len(plist)):
        Pmax[p] = [max(P[(j, p)][kk] for j in range(len(Js))) for kk in range(M + 1)]
    def lhs(Pvec, N):
        return sum(comb(M, kk) * N**(M - kk) * F(Pvec[kk], SCALE) for kk in range(M + 1))
    worst = (F(0), None); refined = 0
    for (m, g, p) in classes:
        N = N_bound(m, g); rhs = Lam * (q**m * Lg[g] * ctab[m][g])**M
        val = lhs(Pmax[p], N)
        if val > rhs:
            refined += 1
            val = max(lhs(P[(j, p)], N) for j in range(len(Js)))
        r = val / rhs
        if r > worst[0]: worst = (r, (m, g))
        assert val <= rhs, f"(C) fails at class m={m} g={g}: ratio {float(r)}"
    # s = 1
    m1, g1 = class_of_digits(dig[T:T+D+200]); N1 = q**T * y[T] * cfun(m1, g1)
    val1 = max(lhs(P[(j, one_pid)], N1) for j in range(len(Js)))
    assert val1 <= Lam, f"(C) fails at s=1: {float(val1/Lam)}"
    print(f"(C) holds for all {len(classes)} classes and s=1; worst ratio {float(worst[0]):.6f} "
          f"at {worst[1]}, s=1 ratio {float(val1/Lam):.6f}; per-J refinement used for {refined} classes")

    # (K)
    assert 3**T * kappa**(M*T) < Lam
    # (Tr): (rB/M) ln 3 <= ((1-B)+(r-1)B) ln(3/2) for r=1..M (implies the condition, as lambda >= 1).
    # Linear in r, so it suffices: r=1, i.e. 3^{B} <= (3/2)^{M(1-B)}, and slope 3^{1/M} <= 3/2.
    b = B.numerator * (1 - B).denominator; c1 = M * (1 - B).numerator * B.denominator
    assert 3**b * 2**c1 <= 3**c1, "(Tr) fails at r=1"
    assert 3 * 2**M <= 3**M, "(Tr) slope fails"
    # exponent: S >= log2(lambda) = log2(Lambda)/T with Lambda <= 2^{T*S}
    if replay:
        S = F(tab["S_upper"])
    else:
        S = F(ceil(float(args.lam_exp) / 0.6931471805599453 * 10**5) + 1, 10**5)
        while not (Lam**S.denominator <= 2**(T * S.numerator)): S += F(1, 10**5)
    assert Lam**S.denominator <= 2**(T * S.numerator)
    assert 3**200 < 2**317                                  # a < 0.585
    assert M >= 2 and F(1, 2) <= B < 1
    assert 3**B.numerator > 2**B.denominator
    saving_factor = 1 - F(585, 1000) - B * S
    theta = (1 - B)**M / (2 - B)
    theta_terminal = (1 - B)**(M - 1) / (2 - B)
    dstar = saving_factor * theta / (1 - theta)
    dstar_terminal = saving_factor * theta_terminal / (1 - theta_terminal)
    cell_mode = "all cells independently" if args.no_symmetry else "reflection enabled"
    print(f"Phase check: {cell_mode}. (K), (Tr) hold. log2(lambda) <= {float(S):.5f}.")
    print(f"Certified delta_* lower bound (rounded), {M} searches: {float(dstar):.7e}")
    print(f"Certified delta_* lower bound (rounded), terminal block ({M - 1} searches): "
          f"{float(dstar_terminal):.7e}")
    if args.save:
        json.dump({"T": T, "M": M, "D": D, "Lambda": str(Lam), "rho": str(rho), "kappa": str(kappa),
                   "B": str(B), "mesh": k, "cut_depth": i0, "cells": len(Js),
                   "cylinders": cyls, "c": [[str(v) for v in row] for row in ctab],
                   "delta_star_lower_M_searches": str(dstar),
                   "delta_star_lower_terminal_block": str(dstar_terminal),
                   "S_upper": str(S)}, open(args.save, "w"))
        print("saved", args.save)

if __name__ == "__main__":
    main()
