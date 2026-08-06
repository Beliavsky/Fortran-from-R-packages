/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   GlmNetCvCpp.hpp
 * Author: Xin Chen
 *
 * Created on July 6, 2017, 5:07 PM
 */

#ifndef GLMNETCVCPP_H
#define GLMNETCVCPP_H
#include <tuple>        // std::make_tuple
#include <algorithm>    // std::shuffle
#include <array>        // std::array
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#include "GlmNetCpp.h"

class GlmNetCvCpp {
public:
        GlmNetCvCpp(const Eigen::MatrixXd& predictor_matrix,
            const Eigen::VectorXd& response_vector,
            double alpha = 1, int num_lambda = 100, int glm_type = 1,
            int max_iter = 1000,
            double abs_tol = 1.0e-4,
            double rel_tol = 1.0e-2,
            bool normalize_grad = false,
            int k_fold = 5,
            bool has_intercept = true,
            int k_fold_iter = 5);
            // function for fitting GLM model given fixed lambda
    Eigen::VectorXd FitGlmFixed();

    // function for generating a grid of candidate lambdas
    Eigen::VectorXd GenerateLambdaGrid();

    // function for automatically choosing the optimal lambda
    // and the corresponding weights using cross validation
    Eigen::VectorXd FitGlmCv();

    // function to compute the smallest lambda that gives zero solution
    double ComputeLambdaMax();

    // Generate training and testing sets for cross validation
    std::tuple<Eigen::MatrixXd, Eigen::MatrixXd,
        Eigen::VectorXd, Eigen::VectorXd> GenerateCvData();
private:
    // predictor_matrix_ is the matrix of the independent variables
    const Eigen::MatrixXd& predictor_matrix_;

    // b is the vector of dependent variables
    const Eigen::VectorXd& response_vector_;

    // alpha is the weight between L1 and L2 regularization, between 0 and 1.
    double alpha_;

    // num_lambda is the number of lambdas for the search grid
    int num_lambda_;

    // type of GLM
    // 1: Exponential
    // 2: Gamma
    int glm_type_;

    // max number of iterations for
    int max_iter_;

    // absolute tolerance
    double abs_tol_;

    // relative tolerance
    double rel_tol_;

    // switch for normalizing the gradient
    bool normalize_grad_;

    // number of folds for cross validation
    int k_fold_;

    // whether the inputs contain intercept or not (as the first variable)
    bool has_intercept_;

    // number of iterations for the cross validation
    int k_fold_iter_;

    // struct for the results of the cross validation for each fixed lambda
//    struct CvResult{
//        double estimated_error;
//        double lambda;
//    };
};

#endif /* GLMNETCVCPP_H */

