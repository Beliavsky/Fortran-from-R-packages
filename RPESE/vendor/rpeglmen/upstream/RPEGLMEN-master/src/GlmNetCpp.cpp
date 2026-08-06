/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "GlmNetCpp.h"   

GlmNetCpp::GlmNetCpp(const Eigen::MatrixXd& predictor_matrix,
        const Eigen::VectorXd& response_vector, double alpha,
        int glm_type,
        int max_iter,
        double abs_tol,
        double rel_tol,
        bool normalize_grad,
        bool has_intercept) :
predictor_matrix_(predictor_matrix),
response_vector_(response_vector),
alpha_(alpha),
glm_type_(glm_type),
max_iter_(max_iter),
abs_tol_(abs_tol),
rel_tol_(rel_tol),
normalize_grad_(normalize_grad),
has_intercept_(has_intercept) {
    //    predictor_matrix_ = A;
    //    response_vector_ = b;
    //    alpha_ = alpha;
    //    num_lambda_ = num_lambda;
    //    glm_type_ = glm_type;
}

double GlmNetCpp::ExpNegativeLogLikelihood(const Eigen::VectorXd& x) {
    // compute the linear component
    Eigen::VectorXd rs = predictor_matrix_ * x;

    // return the negative log-likelihood
    return rs.sum() + response_vector_.transpose() * (-rs).array().exp().matrix();
}

// function to compute the gradient of the negative log-likelihood of exponential GLM

Eigen::VectorXd GlmNetCpp::GradExpNegativeLogLikelihood(const Eigen::VectorXd& x) {
    // number of variables
    // int p = static_cast<int>(predictor_matrix_.cols());

    // number of observations
    int n = static_cast<int>(predictor_matrix_.rows());

    // create vector of n 1s
    Eigen::VectorXd my_ones = Eigen::VectorXd::Ones(n);

    // the gradient of the rs.sum() term in the NLL
    Eigen::VectorXd grad = predictor_matrix_.transpose() * my_ones;

    // compute the linear component
    Eigen::VectorXd rs = predictor_matrix_ * x;

    // the gradient of the response_vector__.transpose() * (-rs).array().exp().matrix() term
    grad += (-predictor_matrix_.transpose()) * ((-rs).array().exp().matrix()).cwiseProduct(response_vector_);

    if (normalize_grad_) {
        return grad.normalized();
    }
    return grad;

}

// function to compute the negative log-likelihood (NLL) of Gamma GLM from data

double GlmNetCpp::GammaNegativeLogLikelihood(const Eigen::VectorXd& x) {
    return 0;
}

// function to compute the gradient of the negative log-likelihood of Gamma GLM

Eigen::VectorXd GlmNetCpp::GradGammaNegativeLogLikelihood(const Eigen::VectorXd& x) {
    return Eigen::VectorXd::Zero(static_cast<int>(predictor_matrix_.cols()));
}

// function for the prox of L1 regularizer,
// which is soft-thresholding operator,
// this is multi-dimensional

Eigen::VectorXd GlmNetCpp::prox_L1(const Eigen::VectorXd& x, double threshold) {
  Eigen::VectorXd tmp = ((abs(x.array()) - threshold).max(0) * x.array().sign()).matrix();
  if (has_intercept_){
    tmp(0) = x(0);
  }
    return tmp;
}

// function for the L1 regularizer
double GlmNetCpp::regularizer_L1(const Eigen::VectorXd& x){
  if(has_intercept_){
    return x.tail(x.size()-1).lpNorm<1>();
  }
    return x.lpNorm<1>();
}

// function for the ENet regularizer
double GlmNetCpp::regularizer_ENet(const Eigen::VectorXd& x){
  if(has_intercept_){
    return alpha_ * regularizer_L1(x) + (1 - alpha_) / 2 * x.tail(x.size()-1).squaredNorm();
  }
    return alpha_ * regularizer_L1(x) + (1 - alpha_) / 2 * x.squaredNorm();
}

// function to compute the value of the entire objective function
// f(x) + lambda * regularizers
double GlmNetCpp::ObjFun(const Eigen::VectorXd& x, double lambda){
    return SmoothObjFun(x) + lambda * regularizer_ENet(x);
}

// function for the smooth part of the objective function

double GlmNetCpp::SmoothObjFun(const Eigen::VectorXd& x) {
    if (glm_type_ == 1)
        return ExpNegativeLogLikelihood(x);
    if (glm_type_ == 2)
        return GammaNegativeLogLikelihood(x);
    return 0;
}

// function for the gradient of the smooth part of the objective function

Eigen::VectorXd GlmNetCpp::GradSmoothObjFun(const Eigen::VectorXd& x) {
    if (glm_type_ == 1)
        return GradExpNegativeLogLikelihood(x);
    if (glm_type_ == 2)
        return GradGammaNegativeLogLikelihood(x);
    return Eigen::VectorXd::Zero(static_cast<int>(predictor_matrix_.cols()));
}

// function for performing Proximal Gradient Descent (PGD)
// This is in fact based on the Fast Proximal gradient at
// https://web.stanford.edu/~boyd/papers/prox_algs/lasso.html#9

Eigen::VectorXd GlmNetCpp::ProxGradDescent(double lambda) {

    double t = 1;
    double beta = 0.5;

    if (rel_tol_ < 0 )
      return Eigen::VectorXd::Zero(static_cast<int>(predictor_matrix_.cols()));

    int num_params = static_cast<int>(predictor_matrix_.cols());

    Eigen::VectorXd x = Eigen::VectorXd::Zero(num_params);
    Eigen::VectorXd xprev = Eigen::VectorXd::Zero(num_params);
    Eigen::VectorXd z;
    int k = 0;
    double obj_val = 0;
    double obj_val_prev = 0;

    while (k < max_iter_) {
//        Eigen::VectorXd y = x + (k / (k + 3)) * (x - xprev);
        Eigen::VectorXd y = x;
//         std::cout << "y = " << y << std::endl;
        while (1) {
            Eigen::VectorXd grad_y = GradSmoothObjFun(y);

//            std::cout << "grad_y =" << grad_y << std::endl;

            double threshold = t * lambda * alpha_
                                / (1 + t * lambda * (1 - alpha_));

            z = prox_L1(y - t * grad_y, threshold);
//            std::cout << "z = " << z << std::endl;

            double lhs = SmoothObjFun(z);
//            std::cout << "lhs =" << lhs << std::endl;

            double rhs1 = SmoothObjFun(y);
//            std::cout << "rhs1 =" << rhs1 << std::endl;

            double rhs2 = grad_y.transpose() * (z - y);
//            std::cout << "rhs2 =" << rhs2 << std::endl;

            double rhs3 = (1 / (2 * t)) * (z - y).squaredNorm();
//            std::cout << "rhs3 =" << rhs3 << std::endl;
            //
            //            std::cout << (lhs <= rhs1 +rhs2 + rhs3) << std::endl;

            if (lhs <=
                    rhs1 +
                    rhs2 +
                    rhs3) {
//                std::cout << "breaking out" << std::endl;
                break;
            }


//            std::cout << "t = " << t << std::endl;
//           std::cout << "beta = " << beta << std::endl;
            t = beta * t;

        }
        xprev = x;
        x = z;

        obj_val_prev = obj_val;
        obj_val = ObjFun(x, lambda);

        if (k > 1)
            if (fabs(obj_val - obj_val_prev) < abs_tol_)
                break;
        k++;
    }
//    std::cout << "num_iter = " << k << std::endl;
    return x;
}

// function to computed the predicted response_vector given the vector of coefficients
// and the predictor_matrix
Eigen::VectorXd GlmNetCpp::Predict(const Eigen::MatrixXd& predictor_matrix_test,
        const Eigen::VectorXd& training_coefficients){
    Eigen::VectorXd res = (predictor_matrix_test * training_coefficients).array().exp();
    return res;
}

// get functions

Eigen::MatrixXd GlmNetCpp::get_predictor_matrix() {
    return predictor_matrix_;
}

Eigen::VectorXd GlmNetCpp::get_response_vector() {
    return response_vector_;
}

double GlmNetCpp::get_alpha() {
    return alpha_;
}

int GlmNetCpp::get_glm_type() {
    return glm_type_;
}

// set functions
//void GlmNetCpp::set_predictor_matrix(Eigen::MatrixXd M){
//    predictor_matrix_ = M;
//}
//
//void GlmNetCpp::set_response_vector(Eigen::VectorXd V){
//    response_vector_ = V;
//}

void GlmNetCpp::set_alpha(double x) {
    alpha_ = x;
}

void GlmNetCpp::set_glm_type(int x) {
    glm_type_ = x;
}


