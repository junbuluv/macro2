function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
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

assert(length(T) >= 11);

T(1) = y(1)^(-params(6));
T(2) = params(1)*params(5)^(1-params(6));
T(3) = 1-(params(3)+exp(y(5))*params(12)*(y(3)-1)+params(13)/2*(y(3)-1)^2);
T(4) = y(3)^params(2);
T(5) = y(4)^(params(2)-1);
T(6) = y(2)^(1-params(2));
T(7) = y(4)^params(2);
T(8) = y(2)^(-params(2));
T(9) = T(4)*exp(y(6))*(1-params(2))*T(7)*T(8);
T(10) = params(2)*exp(y(6))*y(3)^(params(2)-1);
T(11) = T(5)*T(10);

end
