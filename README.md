# The Schneider non-convexity index of Minkowski sums on the line

Yitzchak Shmalo, Einstein Institute of Mathematics, The Hebrew University of Jerusalem.

The Schneider non-convexity index of a compact set A is the least λ ≥ 0 for which A + λ conv A is convex. For compact sets A_1, ..., A_k on the line with prescribed indices c_1 ≤ ... ≤ c_k, the index of A_1 + ... + A_k is at least (1 - B)_+ / (2^k - 1 + B), where B = Σ_i 2^(k-i) (1 - c_i)/(1 + c_i). M. Meyer (Measuring the convexity of compact sumsets with the Schneider non-convexity index, arXiv:2405.00221, Problem 2.2 and Proposition 2.1) constructed sets attaining this value and proved the bound for k = 2 (J. Inequal. Appl. 2025, Paper No. 103). The bound for every k answers his question on the minimal index of a k-fold sum on the line. As a consequence, the least number of compact subsets of the line with common index c < 1 whose sum is convex is 1 + ⌈log2(1/(1 - c))⌉.

## Contents

- `paper/main.tex`, `paper/main.pdf`: the paper.
- `code/check_line_sums.py`: in exact rational arithmetic, checks the inequality on random finite unions of intervals with k = 2, 3, 4 (generic unions, and two-cluster sets whose index is close to 1, for which the bound is mostly positive), checks that Meyer's sets attain it for every sorted index vector on the grid {0, 1/8, ..., 1} with k ≤ 4, and checks 12,000 two-piece families with the same indices at other gap positions and scales. Run `python3 check_line_sums.py`; standard library only, under a minute.
- `check_output.txt`: the output of that run.
- `SHA256SUMS`: sha256 of each file above.
