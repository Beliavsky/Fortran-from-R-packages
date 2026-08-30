#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
build=${1:-"$root/build-strict"}
rm -rf "$build"
mkdir -p "$build/proxy_mod" "$build/proxy_obj" "$build/e1071_mod" "$build/e1071_obj" \
         "$build/classint_mod" "$build/classint_obj" "$build/bin"

flags=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fimplicit-none -fcheck=all \
       -ffpe-trap=invalid,zero,overflow -O0)

proxy_root="$root/dependencies/e1071-fortran/dependencies/proxy-fortran"
proxy_modules=(
    proxy_kinds proxy_ieee proxy_utils proxy_numeric_measures proxy_binary_measures
    proxy_nominal_measures proxy_gower proxy_string_measures proxy_registry proxy_engine
    proxy_dist_utils proxy
)
for module in "${proxy_modules[@]}"; do
    gfortran "${flags[@]}" -J"$build/proxy_mod" -I"$build/proxy_mod" \
        -c "$proxy_root/src/$module.f90" -o "$build/proxy_obj/$module.o"
done

e1071_root="$root/dependencies/e1071-fortran"
e1071_modules=(
    e1071_kinds e1071_constants e1071_array e1071_special e1071_rng e1071_utils
    e1071_discrete e1071_signal e1071_graph e1071_fuzzy e1071_fuzzy_api
    e1071_naive_bayes e1071_knn e1071_lca e1071_ica e1071_matching
    e1071_fclust_indices e1071_svm_solver e1071_svm e1071_sparse e1071_svm_sparse e1071_svm_io e1071_bclust
    e1071_probplot e1071_tune e1071
)
for module in "${e1071_modules[@]}"; do
    gfortran "${flags[@]}" -J"$build/e1071_mod" -I"$build/e1071_mod" -I"$build/proxy_mod" \
        -c "$e1071_root/src/$module.f90" -o "$build/e1071_obj/$module.o"
done

classint_modules=(
    classint_kinds classint_types classint_utils classint_pretty classint_fft classint_dpih
    classint_fisher classint_cluster classint_core classint_metrics classint_api classint
)
for module in "${classint_modules[@]}"; do
    gfortran "${flags[@]}" -J"$build/classint_mod" -I"$build/classint_mod" -I"$build/e1071_mod" -I"$build/proxy_mod" \
        -c "$root/src/$module.f90" -o "$build/classint_obj/$module.o"
done

objects=("$build"/proxy_obj/*.o "$build"/e1071_obj/*.o "$build"/classint_obj/*.o)
for source in "$root"/test/*.f90; do
    name=$(basename "${source%.f90}")
    gfortran "${flags[@]}" -I"$build/classint_mod" -I"$build/e1071_mod" -I"$build/proxy_mod" "$source" \
        "${objects[@]}" -o "$build/bin/$name"
    "$build/bin/$name"
done
for source in "$root"/example/*.f90; do
    name=$(basename "${source%.f90}")
    gfortran "${flags[@]}" -I"$build/classint_mod" -I"$build/e1071_mod" -I"$build/proxy_mod" "$source" \
        "${objects[@]}" -o "$build/bin/$name"
    "$build/bin/$name"
done

(cd "$root" && python tools/check_source_rules.py)
(cd "$e1071_root" && python tools/check_source_rules.py)
(cd "$proxy_root" && python tools/check_source_rules.py)
