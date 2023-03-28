// Samuel Sánchez Gutiérrez
// samuelhdo010524@gmail.com

// [[Rcpp::plugins("cpp11")]]
#include <iostream>
#include <tuple>
#include <time.h>
#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// Objects Conversion
Rcpp::NumericMatrix arma_to_R(arma::mat m){
  
  Rcpp::NumericMatrix r_m(m.n_rows,m.n_cols);
  std::copy(m.begin(),m.end(),r_m.begin());
  
  return(r_m);
}
 
Rcpp::NumericVector arma_to_R(arma::vec v){
  
  Rcpp::NumericVector r_v(v.size());
  std::copy(v.begin(),v.end(),r_v.begin());
  
  return(r_v);
}
 
 
// Expected Value
double expectedValueRFi_j(int ni,int nk, int j){
		return (((double)(ni+nk+1))/(nk+1))*j;
}

// Variance
double varianceRFi_j(int ni,int nk, int j){
		double f1 = ((double)j)/(nk+1);
		double f2 = (((double)(ni+nk+1))/(nk+2))*ni;
		return f1*(1-f1)*f2;
}

// B^{\hat{F}_i}
double BFi(arma::uvec RFi,int ni,int nk){
		double s = 0;
		for(int j = 0; j<nk; j++){
				s += pow(RFi(j) - expectedValueRFi_j(ni,nk,j+1),2)/varianceRFi_j(ni,nk,j+1);
		}
		return s/nk;
}

// Calculate depths
std::tuple<arma::vec,arma::vec> computeDepths(arma::mat X1,
                                              arma::mat X2,
                                              arma::mat Z,
                                              std::string depth,
                                              int N,int n1,
                                              int n2){
  
    // Extract funtion from ddalpha R package
    std::string rdepth = "depth." + depth;
    Rcpp::Environment Rddalpha = Rcpp::Environment::namespace_env("ddalpha");
    Rcpp::Function Rdepth = Rddalpha[rdepth];
    
    // Depths
    arma::vec d1(N,arma::fill::zeros);
    arma::vec d2(N,arma::fill::zeros);
    Rcpp::NumericMatrix rX1 = arma_to_R(X1);
    Rcpp::NumericMatrix rX2 = arma_to_R(X2);
    
    d1.head(n1) = Rcpp::as<arma::vec>(Rdepth(
      Named("x",rX1),
      Named("data",rX1)
    ));
  
    d1.tail(n2) = Rcpp::as<arma::vec>(Rdepth(
      Named("x",rX2),
      Named("data",rX1)
    ));
  
    d2.head(n1) = Rcpp::as<arma::vec>(Rdepth(
      Named("x",rX1),
      Named("data",rX2)
    ));
  
    d2.tail(n2) = Rcpp::as<arma::vec>(Rdepth(
      Named("x",rX2),
      Named("data",rX2)
    ));
  
		return std::make_tuple(d1,d2);
}

// Test Statistic
double Bstat(arma::mat X1,arma::mat X2,arma::mat Z,
             std::string depth,
             int N,int n1,int n2){
  
    arma::vec d1, d2;
    std::tie(d1,d2) = computeDepths(X1,X2,Z,depth,N,n1,n2);

		// Ranks
		arma::uvec RF1 = 1 + arma::sort_index(arma::sort_index(d1));
		arma::uvec RF2 = 1 + arma::sort_index(arma::sort_index(d2));

		// Statistic
		arma::uvec RF1_ord = arma::sort(RF1.tail(n2));
		arma::uvec RF2_ord = arma::sort(RF2.head(n1));
		double BF1 = BFi(RF1_ord,n1,n2);
		double BF2 = BFi(RF2_ord,n2,n1);
		
		return std::max(BF1,BF2);
}

// baraleShirkeTest
// [[Rcpp::export]]
Rcpp::List baraleShirkeTest(Rcpp::NumericMatrix rX1,
                            Rcpp::NumericMatrix rX2,
                            std::string depth,
                            int B,float alpha,
                            bool returnDepths,
                            bool returnSamples){
		
		// Conversión a Armadillo
		arma::mat X1 = Rcpp::as<arma::mat>(rX1);
    arma::mat X2 = Rcpp::as<arma::mat>(rX2);
		
		arma::mat Z = arma::join_vert(X1,X2);
		int n1 = X1.n_rows;
		int n2 = X2.n_rows;
		int N = n1 + n2;
		int n_dim = X1.n_cols;
		
		// ObservedStatistic
		double B0 = Bstat(X1,X2,Z,depth,N,n1,n2);

		// Remuestreo
		double p;
		Rcpp::RObject out_Bsamples;
		if(returnSamples){
		  
		    arma::vec Bsamples(B, arma::fill::none);
		    clock_t t = clock(); // Progress report
		    for(int i=0;i<B;i++){
		        arma::uvec perm = arma::randperm(N);
		        arma::mat bX1 = Z.rows(perm.head(n1));
		        arma::mat bX2 = Z.rows(perm.tail(n2));
		        arma::mat Zb = Z.rows(perm);
		    		Bsamples(i) = Bstat(bX1,bX2,Zb,depth,N,n1,n2);
		    		
		    		if(clock()-t > 1*CLOCKS_PER_SEC){ // Progress report
		    		  t = clock();
		    		  std::cout << std::setprecision(3) << (float)i/B*100 << "% completed." << std::endl;
		    		}
		    }
		    
		    // Compute p-value and return samples
		    arma::uvec comp = (B0 <= Bsamples);
		    p = ((double)sum(comp))/B;
		    
		    out_Bsamples = arma_to_R(Bsamples);
		    
		}else{
		  
		  int upper_Bsamples = 0;
		    clock_t t = clock(); // Progress report
		    for(int i=0;i<B;i++){
		        arma::uvec perm = arma::randperm(N);
		        arma::mat bX1 = Z.rows(perm.head(n1));
		        arma::mat bX2 = Z.rows(perm.tail(n2));
		        arma::mat Zb = Z.rows(perm);
		        if(Bstat(bX1,bX2,Zb,depth,N,n1,n2) > B0) upper_Bsamples++;
		    		
		    		if(clock()-t > 1*CLOCKS_PER_SEC){ // Progress report
		    		  t = clock();
		    		  std::cout << std::setprecision(3) << (float)i/B*100 << "% completed." << std::endl;
		    		}
		    }
		    
		    // Compute p-value and return samples
		    out_Bsamples = R_NilValue;
		    p = (double)upper_Bsamples/B;
		}

		// Conclude
		std::string concl;
		if (p<alpha){
				concl = "Reject";
		}else{
				concl = "Do not reject";
		}
		
		// Return depths
		Rcpp::RObject out_depths;
		if(returnDepths){
		  
        arma::vec d1, d2;
        std::tie(d1,d2) = computeDepths(X1,X2,Z,depth,N,n1,n2);
        
        arma::mat depths = arma::join_horiz(arma::mat(d1),
                                            arma::mat(d2));
		    out_depths = arma_to_R(depths);
		    
		}else{
		    out_depths = R_NilValue;
		}

		// Display Results
		std::string results;
		results += "Barale & Shirke Test:\n";
		results += "Two multivariate samples rank test for\n scale/location equality\n";
		results += "H_0: mu_1 = mu_2 and Sigma_1 = Sigma_2\n";
		results += std::string(37,'-');
		results += '\n';
		results += "Test statistic: ";
		results += std::to_string(B0);
		results += '\n';
		results += "Aproximated p-value: ";
		results += std::to_string(p);
		results += '\n';
		results += std::string(37,'-');
		results += '\n';
		results += "Depth Measure: ";
		results += depth;
		results += '\n';
		results += "Decision: ";
		results += concl;
		results += " the null hypothesis at\n ";
		results += std::to_string(alpha);
		results += " significance level.\n";
		results += "There were used ";
		results += std::to_string(B);
		results += " iterations to\n aproximate the statistic's distribution.";
		
		

		// Return
		return Rcpp::List::create(
						Rcpp::Named("Statistic") = B0,
						Rcpp::Named("NIter")     = B,
						Rcpp::Named("PValue")    = p,
						Rcpp::Named("n1")        = n1,
						Rcpp::Named("n2")        = n2,
						Rcpp::Named("alpha")     = alpha,
						Rcpp::Named("Depth")     = depth,
						Rcpp::Named("DepthVals") = out_depths,
						Rcpp::Named("Samples")   = out_Bsamples,
            Rcpp::Named("Message")   = results
		      );
}
