// This file is part of PIQP-R.
//
// Copyright (c) 2023 piqp-r authors
//
// This source code is licensed under the BSD 2-Clause License found in the
// LICENSE file in the root directory of this source tree.

#ifndef PIQP_TYPES_H
#define PIQP_TYPES_H

#include <RcppEigen.h>

using Vec = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using Mat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using SparseMat = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;

#endif //PIQP_TYPES_H
