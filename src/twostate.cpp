#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
#include <iterator>
#include <algorithm>


double q_12(double x, double a = 1., double b = 0., double c = 2.) {

	double y;
	y = a*(sin(-2*M_PI*(x-b))) + c;
	y = (y >= 0.0)? y: 0.0;
	return y;
}

double q_12_integrated(double x, double a = 1., double b = 0., double c = 2.) {
  
  double y;
  
  y = a*(cos(-2*M_PI*(x-b)))/(2*M_PI) + c*x;
  y = y - a*cos(2*M_PI*b)/(2*M_PI);
  
  return y;
}
// [[Rcpp::export]]
arma::vec get_beta(int Nbasis, double a, double b, double c) {

	arma::vec time_points = arma::linspace<arma::vec>(0.0,1.0, Nbasis+1);
	time_points = time_points.subvec(0,Nbasis) + 1/((float)Nbasis*2.0);

	arma::vec beta0(Nbasis);
	for(int i=0; i<Nbasis; i++) {
		beta0[i] = q_12(time_points[i], a, b, c)*sqrt(1.0/Nbasis);
	}

	return beta0;
	
}

 
arma::mat Qmat(double t, double a = 1., double b = 0., double c = 2.0) {
  arma::mat Q(2,2);
  Q.zeros();

  Q(0,0) = -(q_12(t,a,b,c));
  Q(0,1) = q_12(t,a,b,c); 

  Q(1,0) = c; 
  Q(1,1) = -c; 

  return Q;
}

const double lambda = 500.0;
const Rcpp::NumericVector StateSpace = Rcpp::NumericVector::create(1., 2.);
// const Rcpp::NumericVector initial_prob(5, 0.20000);

Rcpp::NumericVector Jump_times() {
  
	Rcpp::NumericVector A(1000);
	// double lambda = 100;
	A = Rcpp::rexp(1000, lambda);

	Rcpp::NumericVector temp(1000);
	std::partial_sum(A.begin(), A.end(), temp.begin());

	return temp[temp <= 1.0];
}

double Sim_first_State() {
	Rcpp::NumericVector initial_prob(2, 0.50000);
	Rcpp::NumericVector First = Rcpp::RcppArmadillo::sample(StateSpace, 1, true, initial_prob);
	return First[0];
}

Rcpp::NumericVector Sim_next_State(double i, double JumpT, double a = 1., double b = 0., double c = 1.0) {
  arma::mat Pmat = arma::eye(2,2) + (1/lambda)*Qmat(JumpT,a,b,c);
  int j = floor(i);
  return Rcpp::RcppArmadillo::sample(StateSpace, 1, true, Rcpp::wrap(Pmat.row(j-1)));
}

 
Rcpp::List TwoStateSimulation(double a = 1., double b = 0., double c = 2.0) {
  Rcpp::NumericVector TimePoints = Jump_times();
  int N = TimePoints.size();
  double First = Sim_first_State();
  Rcpp::NumericVector States(N+1, First);
  
  if (N >= 1) {
    for(int i = 0; i< N; i ++) {
      States[i+1] = (Sim_next_State(States[i], TimePoints[i],a,b,c))[0];
    }
  }
  Rcpp::List result;
  result["States"] = States;
  result["JumpTimes"] = TimePoints;
  return result;
}


// [[Rcpp::export]]
Rcpp::List simTWO(double a = 1., 
	double b = 0., 
	double c = 2.0){

	Rcpp::List A;
	A = TwoStateSimulation(a, b, c);

	Rcpp::NumericVector States = A["States"];;
	Rcpp::NumericVector JumpTimes = A["JumpTimes"];
	Rcpp::NumericVector transition_times(0);
	Rcpp::NumericVector enter_times(0);
	Rcpp::NumericVector exit_times(0);

	double s_last = States[States.size()-1];
	JumpTimes.push_front(0.0);
	if (JumpTimes[JumpTimes.size() - 1] != 1.00000) {
		JumpTimes.push_back(1.0);
		States.push_back(s_last);
	}

	int N = States.size();
	int M = JumpTimes.size();
	double exit_flag = 0.0;
	if(N != M) {
		exit_flag = 2.0;
	}
	if(N > 0) {
		for (int i = 0; i < N-1; ++i)
		{
			if(States[i] == 1.0 && States[i+1] == 2.)
			{
				transition_times.push_back(JumpTimes[i+1]);
				exit_times.push_back(JumpTimes[i+1]);
			}
		}
		for (int i = 0; i < N-1; ++i)
		{
			if(States[i] == 2.0 && States[i+1] == 1.)
			{
				enter_times.push_back(JumpTimes[i+1]);
			}
		}
		if(States[0] == 1.) enter_times.push_front(0.0);
		if(States[N-1] == 1. && transition_times[transition_times.size() - 1]!= 1.00) exit_times.push_back(1.0);
	}

	Rcpp::List returnList;
	/*
	if(enter_times.size() != exit_times.size()) {
		exit_flag = 1.0;
		returnList["exit_flag"] = exit_flag;
		returnList["States"] = States;
		returnList["JumpTimes"] = JumpTimes;
		returnList["transition_times"] = transition_times;
		returnList["enter_times"] = enter_times;
		returnList["exit_times"] = exit_times;
	} else{
		returnList["exit_flag"] = exit_flag;
		returnList["transition_times"] = transition_times;
		returnList["enter_times"] = enter_times;
		returnList["exit_times"] = exit_times;
	}*/
	returnList["transition_times"] = transition_times;
	returnList["enter_times"] = enter_times;
	returnList["exit_times"] = exit_times;

	return returnList;
}
