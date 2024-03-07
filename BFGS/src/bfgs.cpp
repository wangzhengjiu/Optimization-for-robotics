#include <iostream>
#include <vector>
#include <cmath>
#include <numeric> 

// Rosenbrock 函数
double rosenbrock(const std::vector<double>& x) {
    double sum = 0.0;
    for (int i = 0; i < x.size() - 1; ++i) {
        double term1 = 100 * std::pow(x[i + 1] - std::pow(x[i], 2), 2);
        double term2 = std::pow(1 - x[i], 2);
        sum += term1 + term2;
    }
    return sum;
}

// Rosenbrock 函数的梯度
std::vector<double> rosenbrock_gradient(const std::vector<double>& x) {
    std::vector<double> grad(x.size(), 0.0);
    for (int i = 0; i < x.size() - 1; ++i) {
        grad[i] += -400 * (x[i + 1] - std::pow(x[i], 2)) * x[i] - 2 * (1 - x[i]);
        grad[i + 1] += 200 * (x[i + 1] - std::pow(x[i], 2));
    }
    return grad;
}

// BFGS 算法
std::vector<double> bfgs_method(const std::vector<double>& initial_guess, double tol, int max_iter) {
    int n = initial_guess.size();
    std::vector<double> x = initial_guess;
    std::vector<std::vector<double>> H(n, std::vector<double>(n, 0.0)); // 初始化 Hessian 逆矩阵估计为单位矩阵
    for (int i = 0; i < n; ++i) {
        H[i][i] = 1.0;
    }

    int iter = 0;
    while (iter < max_iter) {
        std::vector<double> grad = rosenbrock_gradient(x);
        double norm_grad = 0.0;
        for (int i = 0; i < n; ++i) {
            norm_grad += std::pow(grad[i], 2);
        }
        norm_grad = std::sqrt(norm_grad);
        if (norm_grad < tol) {
            break;
        }

        // 计算搜索方向 p
        std::vector<double> p(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                p[i] -= H[i][j] * grad[j];
            }
        }

        // 更新 x
        std::vector<double> x_new = x;
        for (int i = 0; i < n; ++i) {
            x_new[i] += p[i];
        }

        // 计算梯度差
        std::vector<double> grad_new = rosenbrock_gradient(x_new);
        std::vector<double> y(n);
        for (int i = 0; i < n; ++i) {
            y[i] = grad_new[i] - grad[i];
        }

        // 更新 Hessian 逆矩阵估计
        std::vector<double> Hy(n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Hy[i] += H[i][j] * y[j];
            }
        }

        double rho = 1.0 / std::inner_product(y.begin(), y.end(), p.begin(), 0.0);
        double alpha = rho * std::inner_product(p.begin(), p.end(), Hy.begin(), 0.0);

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                H[i][j] += rho * (y[i] * p[j] - Hy[i] * Hy[j]);
            }
        }

        x = x_new;
        iter++;
    }

    return x;
}

int main() {
    std::vector<double> initial_guess = {1.0, 1.0}; // 初始猜测
    double tolerance = 1e-6; // 容忍度
    int max_iterations = 1000; // 最大迭代次数

    std::vector<double> solution = bfgs_method(initial_guess, tolerance, max_iterations);

    std::cout << "Solution: ";
    for (double val : solution) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    std::cout << "Value of Rosenbrock function at solution: " << rosenbrock(solution) << std::endl;

    return 0;
}
