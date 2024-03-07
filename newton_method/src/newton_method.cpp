#include <iostream>
#include <vector>
#include <cmath>

double rosenbrock(const std::vector<double>& x) {
    double sum = 0.0;
    for (int i = 0; i < x.size() - 1; ++i) {
        double term1 = 100 * std::pow(x[i + 1] - std::pow(x[i], 2), 2);
        double term2 = std::pow(1 - x[i], 2);
        sum += term1 + term2;
    }
    return sum;
}

std::vector<double> rosenbrock_gradient(const std::vector<double>& x) {
    std::vector<double> grad(x.size(), 0.0);
    for (int i = 0; i < x.size() - 1; ++i) {
        grad[i] += -400 * (x[i + 1] - std::pow(x[i], 2)) * x[i] - 2 * (1 - x[i]);
        grad[i + 1] += 200 * (x[i + 1] - std::pow(x[i], 2));
    }
    return grad;
}

std::vector<std::vector<double>> rosenbrock_hessian(const std::vector<double>& x) {
    int n = x.size();
    std::vector<std::vector<double>> hessian(n, std::vector<double>(n, 0.0));

    for (int i = 0; i < n - 1; ++i) {
        hessian[i][i] += 1200 * std::pow(x[i], 2) - 400 * x[i + 1] + 2;
        hessian[i][i + 1] += -400 * x[i];
        hessian[i + 1][i] += -400 * x[i];
        hessian[i + 1][i + 1] += 200;
    }
    return hessian;
}

std::vector<double> newton_method(const std::vector<double>& initial_guess, double tol, int max_iter) {
    std::vector<double> x = initial_guess;
    int iter = 0;

    while (iter < max_iter) {
        std::vector<double> grad = rosenbrock_gradient(x);
        std::vector<std::vector<double>> hessian = rosenbrock_hessian(x);

        
        std::vector<double> step(x.size(), 0.0);
        for (int i = 0; i < x.size(); ++i) {
            for (int j = 0; j < x.size(); ++j) {
                step[i] -= hessian[i][j] * grad[j];
            }
        }

        // Update x
        for (int i = 0; i < x.size(); ++i) {
            x[i] += step[i];
        }

        // Check for convergence
        double norm_grad = 0.0;
        for (int i = 0; i < x.size(); ++i) {
            norm_grad += std::pow(grad[i], 2);
        }
        norm_grad = std::sqrt(norm_grad);
        if (norm_grad < tol) {
            break;
        }

        iter++;
    }

    return x;
}

int main() {
    std::vector<double> initial_guess = {1.7, 1.6}; // 
    double tolerance = 1e-6;
    int max_iterations = 1000;

    std::vector<double> solution = newton_method(initial_guess, tolerance, max_iterations);

    std::cout << "Solution: ";
    for (double val : solution) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    std::cout << "Value of Rosenbrock function at solution: " << rosenbrock(solution) << std::endl;

    return 0;
}

