#pragma once
#include "armadillo"
class BWM
{
public:
	arma::cube encode(arma::cube& x, arma::mat& wm);
	arma::cube encode(arma::fcube& x, arma::mat& wm);
	arma::cube encode(arma::cube& x, arma::fmat& wm);
	arma::fcube encode(arma::fcube& x, arma::fmat& wm);

};

