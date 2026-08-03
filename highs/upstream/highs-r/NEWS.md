# **R** HIGHS Interface

## 1.14 2026-05-18
- update bundled HiGHS to v1.14.0
- regenerate `R/constants.R` from `HighsOptions.h` for v1.14.0
- strip `CITATION.cff` from CMakeLists install target (new in v1.13.1) via `inst/patch_highs.R`
- HiGHS v1.14.0 requires C++17; bump `inst/build_highs.sh` accordingly

## 1.12 2025-11-23
- Thank you to Benjamin Schwendinger for fixing the
  `HighsHashTree.h:855:32: runtime error: index 1 out of bounds for type 'NodePtr [1]'`
  in **HiGHS**
  https://github.com/ERGO-Code/HiGHS/pull/2608/commits/c9cbe0b2a283ff0f2aa70a769609ca816189ed0b

## 1.10 2025-04-20
- update to HiGHS v.10
- improve the low-level API
- improve README

## 0.1-11 2024-08-03
- update to HiGHS v1.7.2
- improve API
