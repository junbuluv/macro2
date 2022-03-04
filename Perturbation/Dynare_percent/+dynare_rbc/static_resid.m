function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
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
%   residual
%

if T_flag
    T = dynare_rbc.static_resid_tt(T, y, x, params);
end
residual = zeros(6, 1);
lhs = params(5)*T(1);
rhs = T(2)*T(1)*(T(3)+params(2)*exp(y(6))*T(4)*T(5)*T(6));
residual(1) = lhs - rhs;
lhs = params(4)*(1-exp(y(2)))^(-params(11));
rhs = T(1)*T(9);
residual(2) = lhs - rhs;
lhs = exp(y(5))*params(12)+params(13)*exp(y(3)-1);
rhs = T(6)*T(11);
residual(3) = lhs - rhs;
lhs = exp(y(1))+params(5)*exp(y(4));
rhs = T(3)*exp(y(4))+T(6)*T(7)*exp(y(6))*T(4);
residual(4) = lhs - rhs;
lhs = y(6);
rhs = y(6)*params(7)+x(2);
residual(5) = lhs - rhs;
lhs = y(5);
rhs = y(5)*params(8)+x(1);
residual(6) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
