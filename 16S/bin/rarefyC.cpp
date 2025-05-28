#include <random>
#include <Rcpp.h>
using namespace Rcpp;

//Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
// [[Rcpp::export]]
IntegerVector rarefyC(IntegerVector x, int depth) {
      int n = x.size();
      int sum = 0;
      for(int i=0; i < n; i++){sum+=x[i];}
      bool rev_sample = false;
      if (sum < depth) {
            IntegerVector rff(n, NumericVector::get_na());
            return rff;
      } else if (sum == depth) {
            return x;
      } else if (sum/2 < depth) {
            rev_sample = true; 
            depth = sum-depth;
      }
      
      int tags[sum];
      bool is_sample[sum];
      int count = 0;
      for (int i = 0; i < n; i++) {
            for(int j = 0; j < x[i]; j++) { 
                  tags[count] = i;        //将下标重复x[i]次
                  is_sample[count] = false;
                  count++;
            }
      }
      
      std::random_device rd;
      std::mt19937 gen(rd());
      std::uniform_int_distribution<> dis(0, sum-1);
      for (int i = 0; i < depth; i++) {
            int idx = dis(gen);
            if (is_sample[idx]) {
                  i--;
            } else { 
                  is_sample[idx] = true;
            }
            
      }
      
      IntegerVector rff(n, 0);
      for (int i = 0; i < sum; i++) {
            if (is_sample[i]) {rff[tags[i]]++;}
      }
      
      
      if (rev_sample) {
            for (int i = 0; i < n; i++) {rff[i] = x[i]-rff[i];}
      }
      return rff;
}


/*** R

*/
