#include <iostream>
#include <fstream>
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
double AmijioStepSize(const std::vector<double>& x,const std::vector<double>& grad, 
                      const double& aplha = 1.0,const double& beta = 0.5){
   double step_size = 1.0;
   double phi_0 = RosenBrockFunction(x);
   std::vector<double> x_next(x.size());
   while(1){
      for(int i = 0; i < x.size(); ++i){
         x_next[i] = x[i] - step_size * grad[i];
      }
      double phi_next = RosenBrockFunction(x_next); 
      if(phi_next <= phi_0 - aplha * step_size * std::hypot(grad[0],grad[1]) ){
         break;
      }
      step_size *= beta;
   }
   return step_size;
}
std::vector<double> 
     GradientDescent(const double & init_x, const double & init_y, 
                     const double & learning_rate,
                     const double & tol,
                     const int    & max_iter){
    double x = init_x;
    double y = init_y; 
    double step_size  = learning_rate;
    int iter = 0;             
    std::ofstream out_file("linear_search_2d_amijio.txt",std::ios::out);
    if(!out_file){
      std::cout << "Error.Failed to open outFIle" << std::endl;
    }
    while(iter < max_iter){
        std::vector<double> grad = RosenBrockG({x,y});
        double norm = std::hypot(grad[0],grad[1]);
        if(norm < tol){
            break;  
        }
        out_file << "Iteration\t" << iter
                 << "\tx\t"       << x
                 << "\ty\t"       << y
                 << "\tFuntion value" << RosenBrockFunction({x,y})
                 << "\tgradx\t"   << grad[0]
                 << "\tgrady\t"   << grad[1]
                 << std::endl;
        step_size = AmijioStepSize({x,y},grad);

        x -= step_size * grad[0];
        y -= step_size * grad[1];
        iter ++;
    }
    return {x,y};
}


int main(){
   double x = 2.0, y = 4.0;
   double learning_rate = 0.5;
   double tol   = 1e-9;
   int max_iter = 10000;

   std::vector<double> res 
                         =GradientDescent( x, y,
                                           learning_rate,
                                           tol, max_iter);
   std::cout << " RosenBrock value is " << RosenBrockFunction({res[0],res[1]}) 
             << " x is " << res[0] << " y is " << res[1] << std::endl;

//    std::cout << " x Gradient value of RosenBrock is " << << std::endl;
   return 0;
}
