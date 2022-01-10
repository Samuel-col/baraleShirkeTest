// [[Rcpp::plugins("cpp11")]]
#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// Mahalanobis Distance
double MahDist(vec Obs,vec xBar,mat Sinv){
		vec dif = Obs - xBar;
		return as_scalar(dif.t() * Sinv * dif);
}

// Mahalanobis Depth
double mahDepth(vec Obs,vec xBar, mat Sinv){
		double mahDist = MahDist(Obs,xBar,Sinv);
		return 1/(1+mahDist);
}

vec MahDepth(mat Obs,vec xBar,mat Sinv){
		int n = Obs.n_rows;
		vec depths(n,fill::none);
		for(int i=0;i<n;i++)
				depths(i) = mahDepth(Obs.row(i).as_col(),xBar,Sinv);
		return depths;
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
double BFi(uvec RFi,int ni,int nk){
		double s = 0;
		for(int j = 0; j<nk; j++){
				s += pow(RFi(j) - expectedValueRFi_j(ni,nk,j+1),2)/varianceRFi_j(ni,nk,j+1);
		}
		return s/nk;
}

// Test Statistic
double Bstat(mat X1,mat X2,mat Z,int n1,int n2){

		// Shape/Scale parameters
		vec m1 = mean(X1).as_col();
		vec m2 = mean(X2).as_col();

		mat S1 = cov(X1);
		mat S2 = cov(X2);

		// Depths
		vec d1 = MahDepth(Z,m1,S1.i());
		vec d2 = MahDepth(Z,m2,S2.i());
		
		// Ranks
		uvec RF1 = sort_index(sort_index(d1));
		uvec RF2 = sort_index(sort_index(d2));

		// Statistic
		uvec RF1_ord = sort(RF1.tail(n2));
		uvec RF2_ord = sort(RF2.head(n1));
		double BF1 = BFi(RF1_ord,n1,n2);
		double BF2 = BFi(RF2_ord,n2,n1);
		
		return max(BF1,BF2);
}

// baraleShirkeTest
// [[Rcpp::export]]
List baraleShirkeTest(NumericMatrix rX1, NumericMatrix rX2,int B){
		
		// Conversión a Armadillo
		mat X1 = as<mat>(rX1);
		mat X2 = as<mat>(rX2);
		
		if (X1.n_cols!=X2.n_cols){
		  cout << "ERROR: Las poblaciones no tienen la misma cantidad de dimensiones." << endl;
      return List::create();
		}

		mat Z = join_vert(X1,X2);
		int n1 = X1.n_rows;
		int n2 = X2.n_rows;
		int N = n1 + n2;
		
		// ObservedStatistic
		double B0 = Bstat(X1,X2,Z,n1,n2);

		// Bootstrap
		vec Bsamples(B,fill::none);
		for(int i=0;i<B;i++){
				uvec perm = randperm(N);
				mat bX1 = Z.rows(perm.head(n1));
				mat bX2 = Z.rows(perm.tail(n2));
				Bsamples(i) = Bstat(bX1,bX2,Z,n1,n2);
		}

		// pValue
		uvec comp = B0 >= Bsamples;
		double p = ((double)sum(comp))/B;
		
		string concl;
		if (p<0.05){
				concl = "Reject";
		}else{
				concl = "Do not reject";
		}

		// Display Results
		cout << "Barale & Shirke Test:" << endl;
		cout << "Two multivariate samples rank test for scale/location equality" << endl;
		cout << "H_0: mu_1=mu_2 and Sigma_1=Sigma_2" << endl;
		cout << string(75,'-') << endl;
		cout << "Test statistic: " << B0 << endl;
		cout << "Aproximated p-value: "  << p << endl;
		cout << string(75,'-') << endl;
		cout << "Decision: " << concl << " the null hypothesis at 0.05 significance level." << endl;
		cout << "They were used " << B << " iterations to aproximate the statistic's distribution." << endl << endl << endl;

		// Return
		return List::create(
						Named("Statistic") = B0,
						Named("NIter") = B,
						Named("PValue") = p);
}
