#!/usr/bin/env bash
set -euo pipefail
mode=${1:-debug}
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
fc=${FC:-gfortran}
common=(-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace -ffree-line-length-none)
case "$mode" in
  debug) flags=(-O0 -g -fcheck=all) ;;
  release) flags=(-O2) ;;
  *) echo "mode must be debug or release" >&2; exit 2 ;;
esac
b="build/$mode"
rm -rf "$b"
mkdir -p "$b/mod" "$b/obj" "$b/bin"
mods=(
  robustbase_kinds robustbase_sort robustbase_linalg robustbase_probability robustbase_pca
  robustbase_scale robustbase_psi robustbase_medcouple robustbase_covariance robustbase_detmcd
  robustbase_regression robustbase_fastlts robustbase_fast_algorithms robustbase_bylogreg robustbase_lmrob robustbase_glmrob
  robustbase_nonlinear robustbase_nlrob_methods robustbase_adjoutlyingness robustbase_inference robustbase_utilities robustbase_io robustbase
)
for m in "${mods[@]}"; do
  "$fc" "${common[@]}" "${flags[@]}" -J"$b/mod" -I"$b/mod" -c "src/$m.f90" -o "$b/obj/$m.o"
done
ar rcs "$b/librobustbase.a" "$b"/obj/robustbase_*.o "$b/obj/robustbase.o"
"$fc" "${common[@]}" "${flags[@]}" -J"$b/mod" -I"$b/mod" -c test/test_support.f90 -o "$b/obj/test_support.o"
objs=("$b"/obj/robustbase_kinds.o "$b"/obj/robustbase_sort.o "$b"/obj/robustbase_linalg.o "$b"/obj/robustbase_probability.o "$b"/obj/robustbase_pca.o "$b"/obj/robustbase_scale.o "$b"/obj/robustbase_psi.o "$b"/obj/robustbase_medcouple.o "$b"/obj/robustbase_covariance.o "$b"/obj/robustbase_detmcd.o "$b"/obj/robustbase_regression.o "$b"/obj/robustbase_fastlts.o "$b"/obj/robustbase_fast_algorithms.o "$b"/obj/robustbase_bylogreg.o "$b"/obj/robustbase_lmrob.o "$b"/obj/robustbase_glmrob.o "$b"/obj/robustbase_nonlinear.o "$b"/obj/robustbase_nlrob_methods.o "$b"/obj/robustbase_adjoutlyingness.o "$b"/obj/robustbase_inference.o "$b"/obj/robustbase_utilities.o "$b"/obj/robustbase_io.o "$b"/obj/robustbase.o)
for t in test_scales test_covariance test_regression test_next_batch test_remaining; do
  "$fc" "${common[@]}" "${flags[@]}" -J"$b/mod" -I"$b/mod" "test/$t.f90" "$b/obj/test_support.o" "${objs[@]}" -llapack -lblas -o "$b/bin/$t"
done
for a in demo_robustbase fit_csv; do
  "$fc" "${common[@]}" "${flags[@]}" -J"$b/mod" -I"$b/mod" "app/$a.f90" "${objs[@]}" -llapack -lblas -o "$b/bin/$a"
done
"$fc" "${common[@]}" "${flags[@]}" -J"$b/mod" -I"$b/mod" example/nonlinear_example.f90 "${objs[@]}" -llapack -lblas -o "$b/bin/nonlinear_example"
"$fc" "${common[@]}" "${flags[@]}" -J"$b/mod" -I"$b/mod" example/next_batch_example.f90 "${objs[@]}" -llapack -lblas -o "$b/bin/next_batch_example"
