# Certificates for the \(3/2\) counting bounds

This package accompanies *A power saving in the count of generalized \(Z\)-numbers*. It includes the two certificates used in the paper. The verifier programs use Python 3's standard library and a C++20 compiler with 128-bit integer support (for example, Clang or GCC).

| File | Role |
| --- | --- |
| `params.json` | Rational potential for the seven-digit renewal window in Appendix A. |
| `exact/prepare.py` | Generates rational interval enclosures and checks the potential and renewal conditions in Appendix A. |
| `exact/check.cpp` | Checks the phase-window inequalities with outward-rounded integer arithmetic. |
| `exact/check_independent.cpp` | Rechecks all width-seven phase-window inequalities by direct digit-word enumeration and separate interval arithmetic. |
| `independent_interval.hpp` | Provides the independent Horner interval calculation and exact square-root checks used in both cross-checks. |
| `tile/cert_T12_M8.json` | Exact potential and parameters for the twelve-digit, eighth-moment certificate in Appendix B. |
| `tile/prep_tile.py` | Builds the exact input for the tile checker and tests the finite moment inequalities over all interval-state classes. |
| `tile/tilecheck.cpp` | Bounds the tile power sums and counts exit-cell multiplicities. Supports independent computation of all cells. |
| `tile/tile_bound.py` | Checks the rational exponent and prefactor bounds, including the terminal-block placement, for Appendix C. |
| `LICENSE` | MIT License for the verification source code. |
| `SHA256SUMS` | SHA-256 hashes for the other files in the archive. |

## License

The Python and C++ verification source files (including the C++ header) listed above are licensed under the [MIT License](LICENSE).

## Verify the files

Run the commands below from the extracted `cert_age` directory.

~~~sh
shasum -a 256 -c SHA256SUMS
~~~

## Seven-digit window (Appendix A)

In Bash or Zsh:

~~~sh
set -o pipefail
build=$(mktemp -d)
c++ -O3 -std=c++20 exact/check.cpp -o "$build/windowcheck"
python3 exact/prepare.py | "$build/windowcheck"
~~~

Expected final line: `CERTIFIED`. The checker covers 486 patterns and 233,388 intervals. It reports a minimum rational window slack of \(351990680/2^{42}\).

To recompute the phase-window inequalities with separate interval arithmetic and direct enumeration of the digit words:

~~~sh
c++ -O3 -std=c++20 exact/check_independent.cpp -o "$build/window-independent"
python3 exact/prepare.py | "$build/window-independent"
~~~

The rational input is shared with the main check. The independent verifier covers all 233,388 intervals and reports a minimum slack of \(90109619714/2^{50}\).

In `params.json`, `kappa` stores the depth-eight cylinder potentials, `lambda` stores the initial \(v_{\ell_0}\), and `K0` stores \(K_{\rm cap}\). The gap cap \(g_{\max}\) is stored as `G`.

## Twelve-digit moment (Appendix B)

~~~sh
build=$(mktemp -d)
c++ -O3 -std=c++20 tile/tilecheck.cpp -o "$build/tilecheck"
python3 tile/prep_tile.py tile/cert_T12_M8.json --bin "$build/tilecheck" \
    --cache "$build/reflection.json"
python3 tile/prep_tile.py tile/cert_T12_M8.json --bin "$build/tilecheck" \
    --no-symmetry --cache "$build/direct.json"
~~~

The runs agree exactly on 32 exit-cell counts and 57,600 integer power sums. To compare the saved outputs:

~~~sh
python3 - "$build/reflection.json" "$build/direct.json" <<'PY'
import json
from pathlib import Path
import sys

reflection, direct = (json.loads(Path(p).read_text()) for p in sys.argv[1:])
assert direct["key"] == reflection["key"] + "\n--no-symmetry"
assert direct["E"] == reflection["E"]
assert direct["P"] == reflection["P"]
print("Cell counts and integer power sums agree.")
PY
~~~

Both runs check 20,628 base classes and the state \(1\), using 200 upper prefixes and 32 phase cells. The first run uses reflection, and the second computes all 32 cells independently. The maximum exit-cell count is 541,586. For condition (R), use the exact fractions in the JSON. The potential multiplier \(\varrho=554795463/200000000\) is the first multiple of \(10^{-9}\) for which \(\Lambda_{\rm t}\varrho^8\ge541586\). The slack is about \(0.00118603\).

The driver prints separately labeled exponent bounds for the eight-search placement and for the terminal-block placement with seven searches. The JSON stores them in `delta_star_lower_M_searches` and `delta_star_lower_terminal_block`.

In this tile JSON, `T` is \(w_{\rm t}\), `M` is \(M_{\rm t}\), `Lambda` is \(\Lambda_{\rm t}\), `rho` is \(\varrho\), and `kappa` is the age-norm parameter \(\kappa_{\rm age}\).

### Integer arithmetic

The driver generates rational enclosures in Python. The C++ checker uses signed 64-bit fixed-point endpoints at scale \(2^{42}\) and signed 128-bit products and power sums. It checks narrowing, modulus inputs, accumulator additions, and the phase grid before using the results. It stops the power recurrence at the moment order \(M_{\rm t}\). Appendix B gives the product bounds and the rules for assigning phase cut points.

To repeat the calculation with Clang's or GCC's undefined-behavior sanitizer, use a fresh output directory and compile as follows:

~~~sh
build=$(mktemp -d)
c++ -O3 -std=c++20 -fsanitize=undefined -fno-sanitize-recover=undefined \
    tile/tilecheck.cpp -o "$build/tilecheck-sanitized"
python3 tile/prep_tile.py tile/cert_T12_M8.json --bin "$build/tilecheck-sanitized"
~~~

The sanitized run reproduced all 32 exit-cell counts and 57,600 integer power sums. The largest scaled sum took 50 bits at order eight and 62 bits over all orders.

### Independent endpoint arithmetic

This build uses the same tile enumeration and final rational comparisons. The endpoint exponentials use Horner recurrence at scale \(2^{50}\), and the modulus rounding is checked by exact integer squaring. Use a fresh run for this comparison.

~~~sh
build=$(mktemp -d)
c++ -O3 -std=c++20 -DCERT_INDEPENDENT_ARITHMETIC \
    tile/tilecheck.cpp -o "$build/tilecheck-independent"
python3 tile/prep_tile.py tile/cert_T12_M8.json \
    --bin "$build/tilecheck-independent" --no-symmetry
~~~

The independent arithmetic passes all finite tests with a largest ratio of \(0.999513\) to six decimal places. The exit-cell counts and power sums from this build agree exactly between the direct and reflected scans.

## Explicit interval bound (Appendix C)

~~~sh
python3 tile/tile_bound.py tile/cert_T12_M8.json --terminal-block \
    --delta 0.00002222 --constant 112000000000000000000000000000
~~~

This checks the prefactor \(1.12\times10^{29}\) and saving \(2.222\times10^{-5}\) by rational arithmetic. The separate rational comparisons for the ineffective saving \(0.01059\) are printed in Appendix C.
