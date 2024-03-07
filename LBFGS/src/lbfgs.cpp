#include <iostream>
#include <vector>
#include <cmath>
#include <numeric> 

double rosenbrock(const std::vector<double>& x) {
    double sum = 0.0;
    for (size_t i = 0; i < x.size() - 1; ++i) {
        double term1 = 100 * std::pow(x[i + 1] - std::pow(x[i], 2), 2);
        double term2 = std::pow(1 - x[i], 2);
        sum += term1 + term2;
    }
    return sum;
}

std::vector<double> rosenbrock_gradient(const std::vector<double>& x) {
    std::vector<double> g(x.size(), 0.0);
    for (size_t i = 0; i < x.size() - 1; ++i) {
        g[i] += -400 * (x[i + 1] - std::pow(x[i], 2)) * x[i] - 2 * (1 - x[i]);
        g[i + 1] += 200 * (x[i + 1] - std::pow(x[i], 2));
    }
    return g;
}
// L-BFGS Algorithm
std::vector<double> lbfgs(const std::vector<double>& initial_guess, double tol, int max_iter, int m) {
    std::vector<double> x = initial_guess; // Initial guess
    std::vector<std::vector<double>> s_hist; // History of search directions
    std::vector<std::vector<double>> y_hist; // History of gradient differences
    std::vector<double> rho_hist; // History of step size factors

    int iter = 0;
    while (iter < max_iter) {
        // Compute gradient
        std::vector<double> grad = rosenbrock_gradient(x);

        // Update history of search directions
        if (s_hist.size() == m) {
            s_hist.erase(s_hist.begin());
            y_hist.erase(y_hist.begin());
            rho_hist.erase(rho_hist.begin());
        }
        if (!s_hist.empty()) {
            std::vector<double> s = x;
            for (size_t i = 0; i < x.size(); ++i) {
                s[i] -= s_hist.back()[i];
            }
            s_hist.push_back(s);
            y_hist.push_back(grad);
            double rho = 1.0 / std::inner_product(s.begin(), s.end(), y_hist.back().begin(), 0.0);
            rho_hist.push_back(rho);
        }

        // Compute search direction
        std::vector<double> p = grad;
        if (!s_hist.empty()) {
            std::vector<double> q = grad;
            for (int i = s_hist.size() - 1; i >= 0; --i) {
                double alpha = rho_hist[i] * std::inner_product(s_hist[i].begin(), s_hist[i].end(), q.begin(), 0.0);
                for (size_t j = 0; j < x.size(); ++j) {
                    q[j] -= alpha * y_hist[i][j];
                }
            }
            p = q;
        }

        // Update step size
        double step_size = 0.01;

        // Update x
        for (size_t i = 0; i < x.size(); ++i) {
            x[i] -= step_size * p[i];
        }

        // Update iteration counter
        iter++;

        // Check for convergence
        if (std::sqrt(std::inner_product(grad.begin(), grad.end(), grad.begin(), 0.0)) < tol) {
            break;
        }
    }

    return x;
}

int main() {
    const int n = 2; 
    std::vector<double> initial_guess = {1.0, 1.0}; 
    double tolerance = 1e-6; 
    int max_iterations = 1000;
    int m = 5; 

    std::vector<double> solution = lbfgs(initial_guess, tolerance, max_iterations, m);

    std::cout << "Solution:";
    for (double val : solution) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    std::cout << "Value of Rosenbrock function at solution: " << rosenbrock(solution) << std::endl;

    return 0;
}
