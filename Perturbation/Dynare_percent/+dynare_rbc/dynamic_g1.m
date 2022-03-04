function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = dynare_rbc.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(6, 16);
g1(1,4)=params(5)*T(19);
g1(1,10)=(-(T(2)*T(9)*T(20)));
g1(1,11)=(-(T(2)*T(3)*T(26)));
g1(1,12)=(-(T(2)*T(3)*T(32)));
g1(1,7)=(-(T(2)*T(3)*T(38)));
g1(1,13)=(-(T(2)*T(3)*(-(exp(y(13))*params(12)*(exp(y(12))-1)))));
g1(1,14)=(-(T(2)*T(3)*T(7)*T(8)));
g1(2,4)=(-(T(13)*T(19)));
g1(2,5)=params(4)*T(21)-T(1)*T(23);
g1(2,6)=(-(T(1)*T(12)*T(11)*(1-params(2))*exp(y(9))*T(27)));
g1(2,1)=(-(T(1)*T(12)*(1-params(2))*exp(y(9))*T(10)*T(33)));
g1(2,9)=(-(T(1)*T(13)));
g1(3,5)=(-(T(16)*T(24)));
g1(3,6)=params(13)*exp(y(6)-1)-T(17)*T(29);
g1(3,1)=(-(T(17)*T(35)));
g1(3,8)=params(12)*exp(y(8));
g1(3,9)=(-(T(16)*T(17)));
g1(4,4)=exp(y(4));
g1(4,5)=(-(T(11)*exp(y(9))*T(10)*T(24)));
g1(4,6)=(-(exp(y(1))*(-(exp(y(6))*params(12)*exp(y(8))+T(4)*exp(y(6))*2*(exp(y(6))-1)))+T(17)*T(11)*exp(y(9))*T(27)));
g1(4,1)=(-(T(18)+T(17)*exp(y(9))*T(10)*T(33)));
g1(4,7)=params(5)*exp(y(7));
g1(4,8)=(-(exp(y(1))*(-(params(12)*exp(y(8))*(exp(y(6))-1)))));
g1(4,9)=(-(T(17)*T(11)*exp(y(9))*T(10)));
g1(5,3)=(-params(7));
g1(5,9)=1;
g1(5,16)=(-1);
g1(6,2)=(-params(8));
g1(6,8)=1;
g1(6,15)=(-1);

end
