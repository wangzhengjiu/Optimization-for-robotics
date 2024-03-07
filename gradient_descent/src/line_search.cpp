#include <iostream>
#include <cmath>
#include <vector>

double RosenBrockFunction(const std::vector<double>& x){

   double res = 0.0;
   int num = x.size();
   for(int i = 0; i <= num/2 ; i++){
      res += 100.0*(x[2*i - 1] - x[2*i])*(x[2*i - 1] - x[2*i])
             + (x[2*i] - 1)*(x[2*i] - 1);
   }
   return res;
}

double GetLearningRate(){
    double alpha = 0.0;

return alpha;
}


int main(){

   return 0;
}
