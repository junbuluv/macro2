function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_g1_tt(T, y, x, params, steady_state, it_)
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

assert(length(T) >= 38);

T = dynare_rbc.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(19) = exp(y(4))*getPowerDeriv(exp(y(4)),(-params(6)),1);
T(20) = exp(y(10))*getPowerDeriv(exp(y(10)),(-params(6)),1);
T(21) = (-exp(y(5)))*getPowerDeriv(1-exp(y(5)),(-params(11)),1);
T(22) = exp(y(5))*getPowerDeriv(exp(y(5)),(-params(2)),1);
T(23) = (1-params(2))*exp(y(9))*T(10)*T(11)*T(22);
T(24) = exp(y(5))*getPowerDeriv(exp(y(5)),1-params(2),1);
T(25) = exp(y(11))*getPowerDeriv(exp(y(11)),1-params(2),1);
T(26) = T(7)*T(25);
T(27) = exp(y(6))*getPowerDeriv(exp(y(6)),params(2),1);
T(28) = exp(y(6))*getPowerDeriv(exp(y(6)),params(2)-1,1);
T(29) = T(15)*params(2)*exp(y(9))*T(28);
T(30) = exp(y(12))*getPowerDeriv(exp(y(12)),params(2),1);
T(31) = T(6)*params(2)*exp(y(14))*T(30);
T(32) = T(8)*T(31)-(exp(y(13))*params(12)*exp(y(12))+T(4)*exp(y(12))*2*(exp(y(12))-1));
T(33) = exp(y(1))*getPowerDeriv(exp(y(1)),params(2),1);
T(34) = exp(y(1))*getPowerDeriv(exp(y(1)),params(2)-1,1);
T(35) = T(14)*T(34);
T(36) = exp(y(7))*getPowerDeriv(exp(y(7)),params(2)-1,1);
T(37) = T(5)*T(36);
T(38) = T(8)*T(37);

end
