#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

double RosenBrockFunction(const std::vector<double>& x,const double& a = 100,const double& b = 1){
   
   double sum = 0.0;

   double term1, term2;
   
      term1 = a * (x[0] * x[0] - x[1])*(x[0] * x[0] - x[1]);
      term2 = (x[0] - b) *(x[0] - b);
      sum   = term1 + term2;
   
   return sum;

}

std::vector<double> RosenBrockG(const std::vector<double>& x,const double& a = 100,const double& b = 1 ){
   std::vector<double>  res(2,0);
   
   double term1, term2;
   term1 = a * 4 * x[0] *(x[0] * x[0] - x [1])  + 2 * (x[0] - b);
   term2 = -2* a * (x[0] * x[0] - x [1]);
   res[0] = term1;
   res[1] = term2;
   return res;
}

std::vector<double> 
     GradientDescent(const double & init_x, const double & init_y, 
                     const double & learning_rate,
                     const double & tol,
                     const int    & max_iter){
    double x = init_x;
    double y = init_y;   
    int iter = 0;             
    std::ofstream outFile("basic_linear_search.txt",std::ios::out);
    if(!outFile){
      std::cout << "ERROR Failed to Open the file" <<std::endl;
    }
    while(iter < max_iter){
        std::vector<double> grad = RosenBrockG({x,y});
        double norm = std::hypot(grad[0],grad[1]);
      //   std::cout << "norm " << norm << std::endl;
        outFile << "Iteration\t" << iter
                << std::setprecision(3)
                << "\tx\t" << x << "\ty\t" << y
                << "\tFunction value\t" << RosenBrockFunction({x,y})
                << "\tgradx\t" << grad[0]
                << "\tgrady\t" << grad[1]
                << std::endl;
        if(norm < tol){

            break;  
        }
      //   std::vector<double> grad = RosenBrockG({x,y});

        x -= learning_rate * grad[0];
        y -= learning_rate * grad[1];
        iter ++;
    }
    return {x,y};
}


int main(){
   double x = 1.5, y = 1.5;
   double learning_rate = 0.001;//学习率要设置的很小 不然 会震荡的很厉害
   double tol   = 1e-9;
   int max_iter = 10000;
   
   std::vector<double> res 
                         =GradientDescent( x, y,
                                           learning_rate,
                                           tol, max_iter);
   std::cout << "RosenBrock value is " << RosenBrockFunction({res[0],res[1]}) 
             << " x is " << res[0] << " y is " << res[1] << std::endl;

   return 0;
}
