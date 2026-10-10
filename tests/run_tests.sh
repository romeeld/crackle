#!/bin/bash
# Crackle regression tests.
#
# For a set of physical regimes, evolves a gas parcel with one crackle call per system timestep and
# compares T against a converged reference (the same code with each timestep split into 100 calls
# at accuracy=0.01). Fails if the mean |ln(T/T_ref)| over the run exceeds the tolerance for that case.
# Also runs targeted checks of specific fixes, and optionally prints Fortran grackle for comparison.
#
# usage: ./run_tests.sh [-g]     (-g: also show original Fortran grackle results; informational only)
# env:   CRACKLE_DATA_DIR  directory with CloudyData_*.h5 (default: <grackle>/input)
#        GRACKLE_CLIB      grackle src/clib directory with the built library (default: ../..)
#
# Requires grackle to have been built (`make` in src/clib) with the crackle sources.

set -u
cd "$(dirname "$0")"
export CRACKLE_DATA_DIR=${CRACKLE_DATA_DIR:-$(cd ../../../.. && pwd)/input}
SHOW_GRACKLE=0; [ "${1:-}" = "-g" ] && SHOW_GRACKLE=1

make -s crackle_harness || { echo "Failed to build crackle_harness"; exit 1; }
H=./crackle_harness
NFAIL=0
TMP=$(mktemp -d); trap 'rm -rf "$TMP"' EXIT

pass() { printf "  %-6s %s\n" PASS "$1"; }
fail() { printf "  %-6s %s\n" FAIL "$1"; NFAIL=$((NFAIL+1)); }

echo "=== Convergence tests: single step (accuracy=0.1) vs 100x substepped reference (accuracy=0.01)"
echo "    metric: mean|ln(T/Tref)| over all steps, and the worst step"
# name  nH  T0  Z  uvb  isrf  nsteps  dt[Myr]  ssm  h2ss  dustevol  tolerance
CASES=(
  "hotCGM      0.01  1e6  0.1 1 -1  6 10   0 0 0  0.20"
  "hotdense    0.1   1e7  0.3 1 -1  8 5    1 0 0  0.02"
  "selfshield  0.3   1e5  0.3 1 -1  8 5    3 0 0  0.04"
  "UVBheat     0.001 1e3  0.1 1 -1  6 10   0 0 0  0.02"
  "hotISM      1     1e6  1   1  1  8 0.2  3 3 1  0.15"
  "warmISM     1     1e4  1   1  1  8 1    3 3 1  0.05"
  "coldISM     100   100  1   1  1  8 1    3 3 1  0.02"
  "noUVB_ISM   10    5e3  1   0  1  8 1    0 3 1  0.05"
)
for c in "${CASES[@]}"; do
  read name nH T0 Z uvb isrf ns dt ssm h2 de tol <<< "$c"
  args="$nH $T0 $Z $uvb $isrf $ns $dt"
  if ! ACC=0.1 $H 2 $args 1 $ssm $h2 $de > $TMP/run 2> $TMP/err || ! ACC=0.01 $H 2 $args 100 $ssm $h2 $de > $TMP/ref 2>> $TMP/err; then
    fail "$name: harness failed: $(head -1 $TMP/err)"; continue
  fi
  read err maxerr <<< $(paste $TMP/run $TMP/ref | awk 'NR>1{d=log($2/$8); d=d<0?-d:d; s+=d; if(d>m)m=d; n++} END{printf "%.4f %.4f", s/n, m}')
  msg=$(printf "%-11s mean=%.3f (tol %.2f)  worst=%.3f  T_end=%.4g K (ref %.4g K)" $name $err $tol $maxerr $(tail -1 $TMP/run | awk '{print $2}') $(tail -1 $TMP/ref | awk '{print $2}'))
  if awk -v e=$err -v t=$tol 'BEGIN{exit !(e==e && e<t)}'; then pass "$msg"; else fail "$msg"; fi
  if [ $SHOW_GRACKLE = 1 ]; then
    gisrf=$(awk -v x=$isrf 'BEGIN{print (x<0)?0:x}')  # Fortran grackle has no non-ISM (isrf<0) flag
    printf "         grackle T_end=%.4g K\n" $($H 1 $nH $T0 $Z $uvb $gisrf $ns $dt 1 $ssm $h2 $de 2>/dev/null | grep -v "^ *[A-Za-z]" | tail -1 | awk '{print $2}')
  fi
done

echo "=== Targeted checks"
# calculate_temperature before any solve_chemistry call (used to return 1 K)
T=$($H 2 0.01 1e6 0.1 1 -1 0 10 1 0 0 0 2>/dev/null | head -1 | awk '{print $2}')
if awk -v t=$T 'BEGIN{exit !(t>0.99e6 && t<1.01e6)}'; then pass "initial calculate_temperature: T=$T K (expect 1e6)"; else fail "initial calculate_temperature: T=$T K (expect 1e6)"; fi

# H2 formation on dust with dust_chemistry=1 but no dust evolution (used to be frozen at initial xH2=0.1)
x=$($H 2 100 100 1 0 1 8 1 1 0 0 0 2>/dev/null | tail -1 | awk '{print $5}')
if awk -v x=$x 'BEGIN{exit !(x>0.5)}'; then pass "H2 on dust without dust evolution: xH2=$x (expect >0.5)"; else fail "H2 on dust without dust evolution: xH2=$x (expect >0.5)"; fi

# RT H2 photodissociation, with and without UVB (used to be ignored)
for uvb in 0 1; do
  x0=$($H 2 100 100 1 $uvb 1 8 1 1 0 0 1 2>/dev/null | tail -1 | awk '{print $5}')
  x1=$(RTH2=1e-11 $H 2 100 100 1 $uvb 1 8 1 1 0 0 1 2>/dev/null | tail -1 | awk '{print $5}')
  msg="RT H2 dissociation (uvb=$uvb): xH2=$x1 with 1e-11/s vs $x0 without (expect <1e-3 and >0.1)"
  if awk -v a=$x0 -v b=$x1 'BEGIN{exit !(a>0.1 && b<1.e-3)}'; then pass "$msg"; else fail "$msg"; fi
done

echo
if [ $NFAIL -eq 0 ]; then echo "All tests passed."; else echo "$NFAIL test(s) FAILED."; fi
exit $NFAIL
