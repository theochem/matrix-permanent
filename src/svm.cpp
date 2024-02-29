#include <iostream>
#include <string>
#include <vector>
#include <cmath>

#include "svm.h"

// --------------------------------------
// class{HardMargin_SVM} -> constructor
// --------------------------------------
HardMargin_SVM::HardMargin_SVM() {
    // Constructor body is empty
}

// ----------------------------------------
// class{HardMargin_SVM} -> function{dot}
// ----------------------------------------
double HardMargin_SVM::dot(const std::vector<double>& x1, const std::vector<double>& x2) {
    double ans = 0.0;
    if (x1.size() != x2.size()) {
        std::cerr << "Error : Couldn't match the number of elements for inner product." << std::endl;
        std::exit(-1);
    }
    for (size_t i = 0; i < x1.size(); ++i) {
        ans += x1[i] * x2[i];
    }
    return ans;
}

// ------------------------------------------
// class{HardMargin_SVM} -> function{train}
// ------------------------------------------
void HardMargin_SVM::train(const std::vector<std::vector<double>>& class1_data, const std::vector<std::vector<double>>& class2_data, const size_t D, const double lr, const double limit) {
    constexpr double eps = 0.0000001;
    size_t N, Ns;
    double beta = 1.0;

    std::vector<std::vector<double>> x;
    std::vector<int> y;
    std::vector<double> alpha;

    // (1.1) Set class 1 data
    for (const auto& data : class1_data) {
        x.push_back(data);
        y.push_back(1);
    }

    // (1.2) Set class 2 data
    for (const auto& data : class2_data) {
        x.push_back(data);
        y.push_back(-1);
    }

    // (2) Set Lagrange Multiplier and Parameters
    N = x.size();
    alpha = std::vector<double>(N, 0.0);

    // (3) Training
    std::cout << "\n";
    std::cout << "/////////////////////// Training ///////////////////////\n";
    bool judge;
    do {
        judge = false;
        double error = 0.0;

        // (3.1) Update Alpha
        for (size_t i = 0; i < N; ++i) {
            double item1 = 0.0;
            for (size_t j = 0; j < N; ++j) {
                item1 += alpha[j] * y[i] * y[j] * dot(x[i], x[j]);
            }

            double item2 = 0.0;
            for (size_t j = 0; j < N; ++j) {
                item2 += alpha[j] * y[i] * y[j];
            }
            
            double delta = 1.0 - item1 - beta * item2;
            alpha[i] += lr * delta;
            if (alpha[i] < 0.0) {
                alpha[i] = 0.0;
            } else if (std::abs(delta) > limit) {
                judge = true;
                error += std::abs(delta) - limit;
            }
        }

        // (3.2) Update Beta
        double item3 = 0.0;
        for (size_t i = 0; i < N; ++i) {
            item3 += alpha[i] * y[i];
        }
        beta += item3 * item3 / 2.0;

        // (3.3) Output Residual Error
        std::cout << "\rerror: " << error;

    } while (judge);
    std::cout << "\n";
    std::cout << "////////////////////////////////////////////////////////\n";

    // (4.1) Description for support vectors
    Ns = 0;
    xs.clear();
    ys.clear();
    alpha_s.clear();
    for (size_t i = 0; i < N; ++i) {
        if (alpha[i] > eps) {
            xs.push_back(x[i]);
            ys.push_back(y[i]);
            alpha_s.push_back(alpha[i]);
            ++Ns;
        }
    }
    std::cout << "Ns (number of support vectors) = " << Ns << "\n";

    // (4.2) Description for w
    std::cout << "weight = [ ";
    w = std::vector<double>(D, 0.0);
    for (size_t j = 0; j < D; ++j) {
        for (size_t i = 0; i < Ns; ++i) {
            w[j] += alpha_s[i] * ys[i] * xs[i][j];
        }
        std::cout << w[j] << " ";
    }
    std::cout << "]\n";

    // (4.3) Description for b
    b = 0.0;
    for (size_t i = 0; i < Ns; ++i) {
        double sum = 0.0;
        for (size_t j = 0; j < D; ++j) {
            sum += w[j] * xs[i][j];
        }
        b += ys[i] - sum;
    }
    b /= Ns;

    std::cout << "bias = " << b << "\n";
    std::cout << "////////////////////////////////////////////////////////\n\n";
}

// -----------------------------------------
// class{HardMargin_SVM} -> function{test}
// -----------------------------------------
void HardMargin_SVM::test(const std::vector<std::vector<double>>& class1_data, const std::vector<std::vector<double>>& class2_data) {
    correct_c1 = 0;
    for (const auto& data : class1_data) {
        if (g(data) == 1) {
            ++correct_c1;
        }
    }

    correct_c2 = 0;
    for (const auto& data : class2_data) {
        if (g(data) == -1) {
            ++correct_c2;
        }
    }

    accuracy = static_cast<double>(correct_c1 + correct_c2) / static_cast<double>(class1_data.size() + class2_data.size());
}

// --------------------------------------
// class{HardMargin_SVM} -> function{f}
// --------------------------------------
double HardMargin_SVM::f(const std::vector<double>& x) {
    return dot(w, x) + b;    
}

// --------------------------------------
// class{HardMargin_SVM} -> function{g}
// --------------------------------------
double HardMargin_SVM::g(const std::vector<double>& x) {
    double fx = f(x);
    return fx >= 0.0 ? 1 : -1;
}
