"""Checks for the minimal Schneider index of a Minkowski sum of compact sets on the line.

On the line c(A) = (largest gap of A) / (diameter of A); see Lemma 2.1 of the paper. Sets are finite unions of closed
intervals (points allowed) with rational endpoints, and all sums and indices are computed exactly with Fractions.

1. For random families A_1, ..., A_k (k = 2, 3, 4), drawn from two generators (generic unions of intervals, and sets made
   of two clusters at the ends of their hull, whose index is close to 1), check c(A_1 + ... + A_k) >= F(c_1, ..., c_k),
   where F(c) = (1 - B)_+ / (2^k - 1 + B), B = sum_i 2^(k-i) (1 - c_i) / (1 + c_i), c_1 <= ... <= c_k.
2. Check that Meyer's sets A_i = {0, h_i} + [0, h_i r_i], h_i = 2^(k-i), attain F exactly for every sorted index vector
   on the grid {0, 1/8, ..., 1} with k <= 4.
3. For random index vectors, check that two-piece sets with the same indices, other gap positions and other scales never
   give a sum of index below F.

Run: python3 check_line_sums.py   (standard library only, under a minute)
"""
import itertools
import random
from fractions import Fraction as Fr


def norm(iv):
    iv = sorted(iv)
    out = []
    for a, b in iv:
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def idx(A):
    A = norm(A)
    D = A[-1][1] - A[0][0]
    if D == 0:
        return Fr(0)
    g = max([A[i + 1][0] - A[i][1] for i in range(len(A) - 1)] + [Fr(0)])
    return g / D


def msum(A, B):
    return norm([(a + c, b + d) for a, b in A for c, d in B])


def ssum(sets):
    S = sets[0]
    for T in sets[1:]:
        S = msum(S, T)
    return S


def F(cs):
    cs = sorted(cs)
    k = len(cs)
    r = [(1 - c) / (1 + c) for c in cs]
    B = sum(Fr(2) ** (k - i - 1) * r[i] for i in range(k))
    M = 2 ** k - 1
    return max(Fr(0), 1 - B) / (M + B)


def rand_set(rng, den):
    n = rng.randint(1, 4)
    pts = sorted(Fr(rng.randint(0, 12 * den), den) for _ in range(2 * n))
    iv = [(pts[2 * i], pts[2 * i + 1]) for i in range(n)]
    if rng.random() < 0.3:
        iv = [(a, a) if rng.random() < 0.5 else (a, b) for a, b in iv]
    return norm(iv)


def hi_set(rng):
    """Two clusters at the ends of [0, D], each a union of one to three pieces with small holes."""
    D = Fr(rng.randint(8, 40))
    out = []
    for side in (0, 1):
        w = Fr(rng.randint(0, 12), 8)
        base = Fr(0) if side == 0 else D - w
        cuts = sorted(Fr(rng.randint(0, 24), 24) * w for _ in range(2 * rng.randint(0, 2)))
        pts = [Fr(0)] + cuts + [w]
        out += [(base + pts[2 * i], base + pts[2 * i + 1]) for i in range(len(pts) // 2)]
    return norm(out)


def meyer_family(cs):
    cs = sorted(cs)
    k = len(cs)
    out = []
    for i, c in enumerate(cs):
        h = Fr(2) ** (k - i - 1)
        r = (1 - c) / (1 + c)
        out.append(norm([(Fr(0), h * r), (h, h * (1 + r))]))
    return out


def main():
    rng = random.Random(20261126)
    gens = (("generic unions", lambda: rand_set(rng, rng.choice([1, 2, 3]))), ("two-cluster sets", lambda: hi_set(rng)))
    for label, gen in gens:
        for k in (2, 3, 4):
            n = npos = 0
            worst = None
            for _ in range({2: 2000, 3: 1500, 4: 400}[k]):
                sets = [gen() for _ in range(k)]
                cs = [idx(A) for A in sets]
                cS = idx(ssum(sets))
                f = F(cs)
                assert cS >= f, (sets, cs, cS, f)
                n += 1
                if f > 0:
                    npos += 1
                    ratio = cS / f
                    if worst is None or ratio < worst:
                        worst = ratio
            ws = "n/a" if worst is None else f"{float(worst):.4f}"
            print(f"{label}, k={k}: {n} families, c(S) >= F in all; F > 0 in {npos}; min c(S)/F over those = {ws}")
    m = 0
    for k in (1, 2, 3, 4):
        grid = [Fr(j, 8) for j in range(9)]
        for cs in itertools.product(grid, repeat=k):
            if list(cs) != sorted(cs):
                continue
            fam = meyer_family(cs)
            assert [idx(A) for A in fam] == sorted(cs)
            assert idx(ssum(fam)) == F(cs), (cs, idx(ssum(fam)), F(cs))
            m += 1
    print(f"Meyer's sets attain F exactly on {m} sorted index vectors (k <= 4, grid 1/8)")
    t = 0
    for _ in range(400):
        k = rng.choice([2, 3])
        cs = sorted(Fr(rng.randint(3, 8), 8) for _ in range(k))
        fam = meyer_family(cs)
        for _ in range(30):
            new = []
            for c, A in zip(cs, fam):
                D = A[-1][1] - A[0][0]
                g = c * D
                left = Fr(rng.randint(0, 16), 16) * (D - g)
                B = [(Fr(0), left), (left + g, D)]
                scale = Fr(rng.randint(1, 8), 4)
                new.append(norm([(a * scale, b * scale) for a, b in B]))
            assert [idx(A) for A in new] == cs
            assert idx(ssum(new)) >= F(cs)
            t += 1
    print(f"{t} two-piece families with prescribed indices, other gap positions and scales: none below F")


if __name__ == "__main__":
    main()
