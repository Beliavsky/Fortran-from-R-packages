/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   GlmNetCvCpp.cpp
 * Author: Xin Chen
 *
 * Created on July 6, 2017, 5:07 PM
 */

#include "GlmNetCvCpp.h"

GlmNetCvCpp::GlmNetCvCpp(const Eigen::MatrixXd& predictor_matrix,
        const Eigen::VectorXd& response_vector, double alpha,
        int num_lambda, int glm_type,
        int max_iter,
        double abs_tol,
        double rel_tol,
        bool normalize_grad,
        int k_fold,
        bool has_intercept,
        int k_fold_iter) :
predictor_matrix_(predictor_matrix),
response_vector_(response_vector),
alpha_(alpha),
num_lambda_(num_lambda),
glm_type_(glm_type),
max_iter_(max_iter),
abs_tol_(abs_tol),
rel_tol_(rel_tol),
normalize_grad_(normalize_grad),
k_fold_(k_fold),
has_intercept_(has_intercept),
k_fold_iter_(k_fold_iter) {
    //    predictor_matrix_ = A;
    //    response_vector_ = b;
    //    alpha_ = alpha;
    //    num_lambda_ = num_lambda;
    //    glm_type_ = glm_type;
}

// function for fitting GLM model given fixed lambda

Eigen::VectorXd GlmNetCvCpp::FitGlmFixed() {
    return Eigen::VectorXd::Zero(3);
}

// function for generating a grid of candidate lambdas

Eigen::VectorXd GlmNetCvCpp::GenerateLambdaGrid() {
    double lambda_max = GlmNetCvCpp::ComputeLambdaMax();
    return (Eigen::VectorXd::LinSpaced(num_lambda_, log(0.001), log(lambda_max))).array().exp();
}

// function for automatically choosing the optimal lambda
// and the corresponding weights using cross validation

Eigen::VectorXd GlmNetCvCpp::FitGlmCv() {


    //  generate lambda grid
    Eigen::VectorXd lambda_grid = GenerateLambdaGrid();

    // initialize vectors to store cv results
    std::vector<double> predicted_errors;

    //  for each lambda in lambda grid
    for (int i = 0; i < lambda_grid.size(); i++) {

        double error = 0;

        // run k_fold cv
        for (int cv_iter = 0; cv_iter < k_fold_iter_; cv_iter++) {
            //  generate training and testing data sets
            Eigen::MatrixXd predictor_matrix_train;
            Eigen::VectorXd response_vector_train;
            Eigen::MatrixXd predictor_matrix_test;
            Eigen::VectorXd response_vector_test;

            std::tie(predictor_matrix_train,
                    predictor_matrix_test,
                    response_vector_train,
                    response_vector_test) = GenerateCvData();

            // construct an GlmNetCpp object using training data
            GlmNetCpp my_glm(predictor_matrix_train,
                    response_vector_train,
                    alpha_,
                    glm_type_,
                    max_iter_,
                    abs_tol_,
                    rel_tol_,
                    normalize_grad_,
                    has_intercept_);

            // find the optimal coefficients using training data
            Eigen::VectorXd training_coeffs = my_glm.ProxGradDescent(lambda_grid(i));

            // use training_coeffs on the predictor_matrix_test
            // to get the predicted responses
            Eigen::VectorXd response_vector_predicted = my_glm.Predict(predictor_matrix_test,
                    training_coeffs);

            // compute the error of the predicted response vector versus the test response vector
            error += (response_vector_predicted - response_vector_test).norm();

        }

        error = error / k_fold_iter_;

        // save the results
        predicted_errors.push_back(error);
    }

    // find the lambda corresponding to the smallest predicted_error
    int min_pos = 0;
    for(int i = 1; i < predicted_errors.size(); i++){
        if (predicted_errors[min_pos] > predicted_errors[i]){
            min_pos = i;
        }
    }

    double best_lambda = lambda_grid(min_pos);

    // train the model using the best_lambda and entire training set
                GlmNetCpp my_glm(predictor_matrix_,
                    response_vector_,
                    alpha_,
                    glm_type_,
                    max_iter_,
                    abs_tol_,
                    rel_tol_,
                    normalize_grad_,
                    has_intercept_);
                Eigen::VectorXd best_coeffs = my_glm.ProxGradDescent(best_lambda);

    return best_coeffs;
}

// function to compute the smallest lambda that gives zero solution

double GlmNetCvCpp::ComputeLambdaMax() {
    return (((response_vector_.array() - 1).matrix().transpose() * predictor_matrix_).cwiseAbs() / alpha_).maxCoeff();
}

// Generate training and testing sets for cross validation

std::tuple<Eigen::MatrixXd, Eigen::MatrixXd,
Eigen::VectorXd, Eigen::VectorXd> GlmNetCvCpp::GenerateCvData() {
    // number of observations
    int num_obs = static_cast<int>(predictor_matrix_.rows());

    // number of variables
    int num_vars = static_cast<int>(predictor_matrix_.cols());

    // genereate eigen::vector of 0 to num_obs - 1
    Eigen::VectorXi idx1 = Eigen::VectorXi::LinSpaced(num_obs, 0, num_obs - 1);

    // convert eigen::vector to std::vector
    std::vector<int> idx(idx1.data(), idx1.data() + idx1.size());

    // obtain a time-based seed:
    unsigned seed = static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count());

    // shuffle the std::vector
    shuffle(idx.begin(), idx.end(), std::default_random_engine(seed));

    // The training set contains (k_fold_-1)/k_fold_ of the data
    int train_size = (k_fold_ - 1) / k_fold_ * num_obs;
    std::vector<int> train_idx;
    train_idx.resize(train_size);
    for (int i = 0; i < train_size; i++) {
        train_idx[i] = idx[i];
    }

    // Construct the training set of the predictor_matrix_ and response_vector_
    Eigen::MatrixXd predictor_matrix_train;
    Eigen::VectorXd response_vector_train;
    predictor_matrix_train.resize(train_size, num_vars);
    response_vector_train.resize(train_size);
    for (int i = 0; i < train_size; i++) {
        predictor_matrix_train.row(i) = predictor_matrix_.row(train_idx[i]);
        response_vector_train(i) = response_vector_(train_idx[i]);
    }

    // The test set contains 1/k of the data
    int test_size = num_obs - train_size;
    std::vector<int> test_idx;
    test_idx.resize(test_size);
    for (int i = train_size; i < num_obs; i++) {
        test_idx[i] = idx[i];
    }

    // Construct the test set of the predictor_matrix_ and response_vector_
    Eigen::MatrixXd predictor_matrix_test;
    Eigen::VectorXd response_vector_test;
    predictor_matrix_test.resize(test_size, num_vars);
    response_vector_test.resize(test_size);
    for (int i = 0; i < test_size; i++) {
        predictor_matrix_test.row(i) = predictor_matrix_.row(test_idx[i]);
        response_vector_test(i) = response_vector_(test_idx[i]);
    }
    return std::make_tuple(predictor_matrix_train, predictor_matrix_test,
            response_vector_train, response_vector_test);
}


