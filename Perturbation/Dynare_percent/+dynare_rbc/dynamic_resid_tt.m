function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 18);

T(1) = exp(y(4))^(-params(6));
T(2) = params(1)*params(5)^(1-params(6));
T(3) = exp(y(10))^(-params(6));
T(4) = params(13)/2;
T(5) = params(2)*exp(y(14))*exp(y(12))^params(2);
T(6) = exp(y(7))^(params(2)-1);
T(7) = T(5)*T(6);
T(8) = exp(y(11))^(1-params(2));
T(9) = 1-(params(3)+exp(y(13))*params(12)*(exp(y(12))-1)+T(4)*(exp(y(12))-1)^2)+T(7)*T(8);
T(10) = exp(y(6))^params(2);
T(11) = exp(y(1))^params(2);
T(12) = exp(y(5))^(-params(2));
T(13) = (1-params(2))*exp(y(9))*T(10)*T(11)*T(12);
T(14) = params(2)*exp(y(9))*exp(y(6))^(params(2)-1);
T(15) = exp(y(1))^(params(2)-1);
T(16) = T(14)*T(15);
T(17) = exp(y(5))^(1-params(2));
T(18) = exp(y(1))*(1-(params(3)+params(12)*exp(y(8))*(exp(y(6))-1)+T(4)*(exp(y(6))-1)^2));

end
