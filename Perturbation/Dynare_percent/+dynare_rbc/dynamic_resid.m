function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = dynare_rbc.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(6, 1);
lhs = params(5)*T(1);
rhs = T(2)*T(3)*T(9);
residual(1) = lhs - rhs;
lhs = params(4)*(1-exp(y(5)))^(-params(11));
rhs = T(1)*T(13);
residual(2) = lhs - rhs;
lhs = params(12)*exp(y(8))+params(13)*exp(y(6)-1);
rhs = T(16)*T(17);
residual(3) = lhs - rhs;
lhs = exp(y(4))+params(5)*exp(y(7));
rhs = T(18)+T(17)*T(11)*exp(y(9))*T(10);
residual(4) = lhs - rhs;
lhs = y(9);
rhs = params(7)*y(3)+x(it_, 2);
residual(5) = lhs - rhs;
lhs = y(8);
rhs = params(8)*y(2)+x(it_, 1);
residual(6) = lhs - rhs;

end