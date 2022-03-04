function T = static_g1_tt(T, y, x, params)
% function T = static_g1_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 17);

T = dynare_rbc.static_resid_tt(T, y, x, params);

T(12) = exp(y(1))*getPowerDeriv(exp(y(1)),(-params(6)),1);
T(13) = exp(y(2))*getPowerDeriv(exp(y(2)),1-params(2),1);
T(14) = exp(y(5))*params(12)*exp(y(3))+params(13)/2*exp(y(3))*2*(exp(y(3))-1);
T(15) = exp(y(3))*getPowerDeriv(exp(y(3)),params(2),1);
T(16) = exp(y(4))*getPowerDeriv(exp(y(4)),params(2)-1,1);
T(17) = exp(y(4))*getPowerDeriv(exp(y(4)),params(2),1);

end
