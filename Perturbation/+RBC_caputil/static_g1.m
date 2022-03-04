function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   g1
%

if T_flag
    T = RBC_caputil.static_g1_tt(T, y, x, params);
end
g1 = zeros(6, 6);
g1(1,1)=params(5)*T(12)-(T(3)+params(2)*exp(y(6))*T(4)*T(5)*T(6))*T(2)*T(12);
g1(1,2)=(-(T(1)*T(2)*params(2)*exp(y(6))*T(4)*T(5)*T(13)));
g1(1,3)=(-(T(1)*T(2)*(T(6)*T(5)*params(2)*exp(y(6))*T(14)-(exp(y(5))*params(12)+params(13)/2*2*(y(3)-1)))));
g1(1,4)=(-(T(1)*T(2)*T(6)*params(2)*exp(y(6))*T(4)*T(15)));
g1(1,5)=(-(T(1)*T(2)*(-(exp(y(5))*params(12)*(y(3)-1)))));
g1(1,6)=(-(T(1)*T(2)*params(2)*exp(y(6))*T(4)*T(5)*T(6)));
g1(2,1)=(-(T(9)*T(12)));
g1(2,2)=params(4)*(-(getPowerDeriv(1-y(2),(-params(11)),1)))-T(1)*T(4)*exp(y(6))*(1-params(2))*T(7)*getPowerDeriv(y(2),(-params(2)),1);
g1(2,3)=(-(T(1)*T(8)*T(7)*exp(y(6))*(1-params(2))*T(14)));
g1(2,4)=(-(T(1)*T(8)*T(4)*exp(y(6))*(1-params(2))*T(16)));
g1(2,6)=(-(T(1)*T(9)));
g1(3,2)=(-(T(11)*T(13)));
g1(3,3)=params(13)-T(6)*T(5)*params(2)*exp(y(6))*getPowerDeriv(y(3),params(2)-1,1);
g1(3,4)=(-(T(6)*T(10)*T(15)));
g1(3,5)=exp(y(5))*params(12);
g1(3,6)=(-(T(6)*T(11)));
g1(4,1)=1;
g1(4,2)=(-(T(7)*exp(y(6))*T(4)*T(13)));
g1(4,3)=(-(y(4)*(-(exp(y(5))*params(12)+params(13)/2*2*(y(3)-1)))+T(6)*T(7)*exp(y(6))*T(14)));
g1(4,4)=params(5)-(T(3)+T(6)*exp(y(6))*T(4)*T(16));
g1(4,5)=(-(y(4)*(-(exp(y(5))*params(12)*(y(3)-1)))));
g1(4,6)=(-(T(6)*T(7)*exp(y(6))*T(4)));
g1(5,6)=1-params(7);
g1(6,5)=1-params(8);
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
