#include <iostream>
#include <cmath>
#include <vector>

double RosenBrockFunction(const std::vector<double>& x,const double& a = 100,const double& b = 1){
   
   double sum = 0.0;

   double term1, term2;
   
      term1 = a * (x[0] * x[0] - x[1])*(x[0] * x[0] - x[1]);
      term2 = (x[0] - b) *(x[0] - b);
      sum   = term1 + term2;
   
   return sum;

}

std::vector<double> RosenBrockG(const std::vector<double>& x,const double& a = 100,const double& b = 1 ){
   std::vector<double>  res(2);
   int n = x.size() / 2 - 1;
   double term1, term2;
   term1 = a * 4 * x[0] *(x[0] * x[0] - x [1])  + 2 * (x[0] - b);
   term2 = -2* a * (x[0] * x[0] - x [1]);
   res[0] = term1;
   res[1] = term2;
   return res;
}

double RobbinsMonroRule(const int k, const double a, const double b) {
    return a / (b + k);
}
std::vector<double> 
     GradientDescent(const double & init_x, const double & init_y, 
                     const double & learning_rate,
                     const double & tol,
                     const int    & max_iter){
    double x = init_x;
    double y = init_y;   
    double step_size = learning_rate;
    int iter = 0;             
    while(iter < max_iter){
        
        std::vector<double> grad = RosenBrockG({x,y});
        double norm = std::hypot(grad[0],grad[1]);
        if(norm < tol){

            break;  
        }
        step_size = RobbinsMonroRule(iter,1.0,1.0);
        x -= step_size * grad[0];
        y -= step_size * grad[1];
        iter ++;
    }
    return {x,y};
}


int main(){
   double x = 10.0, y = 9.0;
   double learning_rate = 0.1;
   double tol   = 1e-9;
   int max_iter = 10000;

   std::vector<double> res 
                         =GradientDescent( x, y,
                                           learning_rate,
                                           tol, max_iter);
   std::cout << " RosenBrock value is " << RosenBrockFunction({res[0],res[1]}) 
             << " x is " << res[0] << " y is " << res[1] << std::endl;

   return 0;
}
