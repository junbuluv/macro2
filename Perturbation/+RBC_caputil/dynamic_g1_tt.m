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

assert(length(T) >= 37);

T = RBC_caputil.dynamic_resid_tt(T, y, x, params, steady_state, it_);

T(19) = getPowerDeriv(y(4),(-params(6)),1);
T(20) = T(2)*getPowerDeriv(y(10),(-params(6)),1);
T(21) = getPowerDeriv(y(5),(-params(2)),1);
T(22) = (1-params(2))*exp(y(9))*T(10)*T(11)*T(21);
T(23) = getPowerDeriv(y(5),1-params(2),1);
T(24) = getPowerDeriv(y(11),1-params(2),1);
T(25) = T(7)*T(24);
T(26) = getPowerDeriv(y(6),params(2),1);
T(27) = params(2)*exp(y(9))*getPowerDeriv(y(6),params(2)-1,1);
T(28) = T(15)*T(27);
T(29) = params(2)*exp(y(14))*getPowerDeriv(y(12),params(2),1);
T(30) = T(6)*T(29);
T(31) = T(8)*T(30)-(exp(y(13))*params(12)+T(4)*2*(y(12)-1));
T(32) = getPowerDeriv(y(1),params(2),1);
T(33) = getPowerDeriv(y(1),params(2)-1,1);
T(34) = T(14)*T(33);
T(35) = getPowerDeriv(y(7),params(2)-1,1);
T(36) = T(5)*T(35);
T(37) = T(8)*T(36);

end
