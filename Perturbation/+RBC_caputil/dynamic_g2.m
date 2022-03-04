function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
% function g2 = dynamic_g2(T, y, x, params, steady_state, it_, T_flag)
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
%   g2
%

if T_flag
    T = RBC_caputil.dynamic_g2_tt(T, y, x, params, steady_state, it_);
end
g2_i = zeros(94,1);
g2_j = zeros(94,1);
g2_v = zeros(94,1);

g2_i(1)=1;
g2_i(2)=1;
g2_i(3)=1;
g2_i(4)=1;
g2_i(5)=1;
g2_i(6)=1;
g2_i(7)=1;
g2_i(8)=1;
g2_i(9)=1;
g2_i(10)=1;
g2_i(11)=1;
g2_i(12)=1;
g2_i(13)=1;
g2_i(14)=1;
g2_i(15)=1;
g2_i(16)=1;
g2_i(17)=1;
g2_i(18)=1;
g2_i(19)=1;
g2_i(20)=1;
g2_i(21)=1;
g2_i(22)=1;
g2_i(23)=1;
g2_i(24)=1;
g2_i(25)=1;
g2_i(26)=1;
g2_i(27)=1;
g2_i(28)=1;
g2_i(29)=1;
g2_i(30)=1;
g2_i(31)=1;
g2_i(32)=2;
g2_i(33)=2;
g2_i(34)=2;
g2_i(35)=2;
g2_i(36)=2;
g2_i(37)=2;
g2_i(38)=2;
g2_i(39)=2;
g2_i(40)=2;
g2_i(41)=2;
g2_i(42)=2;
g2_i(43)=2;
g2_i(44)=2;
g2_i(45)=2;
g2_i(46)=2;
g2_i(47)=2;
g2_i(48)=2;
g2_i(49)=2;
g2_i(50)=2;
g2_i(51)=2;
g2_i(52)=2;
g2_i(53)=2;
g2_i(54)=2;
g2_i(55)=2;
g2_i(56)=2;
g2_i(57)=3;
g2_i(58)=3;
g2_i(59)=3;
g2_i(60)=3;
g2_i(61)=3;
g2_i(62)=3;
g2_i(63)=3;
g2_i(64)=3;
g2_i(65)=3;
g2_i(66)=3;
g2_i(67)=3;
g2_i(68)=3;
g2_i(69)=3;
g2_i(70)=3;
g2_i(71)=3;
g2_i(72)=3;
g2_i(73)=3;
g2_i(74)=4;
g2_i(75)=4;
g2_i(76)=4;
g2_i(77)=4;
g2_i(78)=4;
g2_i(79)=4;
g2_i(80)=4;
g2_i(81)=4;
g2_i(82)=4;
g2_i(83)=4;
g2_i(84)=4;
g2_i(85)=4;
g2_i(86)=4;
g2_i(87)=4;
g2_i(88)=4;
g2_i(89)=4;
g2_i(90)=4;
g2_i(91)=4;
g2_i(92)=4;
g2_i(93)=4;
g2_i(94)=4;
g2_j(1)=52;
g2_j(2)=154;
g2_j(3)=155;
g2_j(4)=170;
g2_j(5)=156;
g2_j(6)=186;
g2_j(7)=151;
g2_j(8)=106;
g2_j(9)=157;
g2_j(10)=202;
g2_j(11)=158;
g2_j(12)=218;
g2_j(13)=171;
g2_j(14)=172;
g2_j(15)=187;
g2_j(16)=167;
g2_j(17)=107;
g2_j(18)=174;
g2_j(19)=219;
g2_j(20)=188;
g2_j(21)=183;
g2_j(22)=108;
g2_j(23)=189;
g2_j(24)=204;
g2_j(25)=190;
g2_j(26)=220;
g2_j(27)=103;
g2_j(28)=110;
g2_j(29)=215;
g2_j(30)=205;
g2_j(31)=222;
g2_j(32)=52;
g2_j(33)=53;
g2_j(34)=68;
g2_j(35)=54;
g2_j(36)=84;
g2_j(37)=49;
g2_j(38)=4;
g2_j(39)=57;
g2_j(40)=132;
g2_j(41)=69;
g2_j(42)=70;
g2_j(43)=85;
g2_j(44)=65;
g2_j(45)=5;
g2_j(46)=73;
g2_j(47)=133;
g2_j(48)=86;
g2_j(49)=81;
g2_j(50)=6;
g2_j(51)=89;
g2_j(52)=134;
g2_j(53)=1;
g2_j(54)=9;
g2_j(55)=129;
g2_j(56)=137;
g2_j(57)=69;
g2_j(58)=70;
g2_j(59)=85;
g2_j(60)=65;
g2_j(61)=5;
g2_j(62)=73;
g2_j(63)=133;
g2_j(64)=86;
g2_j(65)=81;
g2_j(66)=6;
g2_j(67)=89;
g2_j(68)=134;
g2_j(69)=1;
g2_j(70)=9;
g2_j(71)=129;
g2_j(72)=120;
g2_j(73)=137;
g2_j(74)=69;
g2_j(75)=70;
g2_j(76)=85;
g2_j(77)=65;
g2_j(78)=5;
g2_j(79)=73;
g2_j(80)=133;
g2_j(81)=86;
g2_j(82)=81;
g2_j(83)=6;
g2_j(84)=88;
g2_j(85)=118;
g2_j(86)=89;
g2_j(87)=134;
g2_j(88)=1;
g2_j(89)=8;
g2_j(90)=113;
g2_j(91)=9;
g2_j(92)=129;
g2_j(93)=120;
g2_j(94)=137;
g2_v(1)=params(5)*T(38);
g2_v(2)=(-(T(9)*T(2)*getPowerDeriv(y(10),(-params(6)),2)));
g2_v(3)=(-(T(20)*T(25)));
g2_v(4)=g2_v(3);
g2_v(5)=(-(T(20)*T(31)));
g2_v(6)=g2_v(5);
g2_v(7)=(-(T(20)*T(37)));
g2_v(8)=g2_v(7);
g2_v(9)=(-(T(20)*(-(exp(y(13))*params(12)*(y(12)-1)))));
g2_v(10)=g2_v(9);
g2_v(11)=(-(T(7)*T(8)*T(20)));
g2_v(12)=g2_v(11);
g2_v(13)=(-(T(3)*T(7)*getPowerDeriv(y(11),1-params(2),2)));
g2_v(14)=(-(T(3)*T(24)*T(30)));
g2_v(15)=g2_v(14);
g2_v(16)=(-(T(3)*T(24)*T(36)));
g2_v(17)=g2_v(16);
g2_v(18)=(-(T(3)*T(25)));
g2_v(19)=g2_v(18);
g2_v(20)=(-(T(3)*(T(8)*T(6)*params(2)*exp(y(14))*getPowerDeriv(y(12),params(2),2)-params(13))));
g2_v(21)=(-(T(3)*T(8)*T(29)*T(35)));
g2_v(22)=g2_v(21);
g2_v(23)=(-(T(3)*(-(exp(y(13))*params(12)))));
g2_v(24)=g2_v(23);
g2_v(25)=(-(T(3)*T(8)*T(30)));
g2_v(26)=g2_v(25);
g2_v(27)=(-(T(3)*T(8)*T(5)*getPowerDeriv(y(7),params(2)-1,2)));
g2_v(28)=(-(T(3)*T(37)));
g2_v(29)=g2_v(28);
g2_v(30)=(-(T(3)*(-(exp(y(13))*params(12)*(y(12)-1)))));
g2_v(31)=(-(T(3)*T(7)*T(8)));
g2_v(32)=(-(T(13)*T(38)));
g2_v(33)=(-(T(19)*T(22)));
g2_v(34)=g2_v(33);
g2_v(35)=(-(T(19)*T(12)*T(11)*(1-params(2))*exp(y(9))*T(26)));
g2_v(36)=g2_v(35);
g2_v(37)=(-(T(19)*T(12)*(1-params(2))*exp(y(9))*T(10)*T(32)));
g2_v(38)=g2_v(37);
g2_v(39)=(-(T(13)*T(19)));
g2_v(40)=g2_v(39);
g2_v(41)=params(4)*getPowerDeriv(1-y(5),(-params(11)),2)-T(1)*(1-params(2))*exp(y(9))*T(10)*T(11)*getPowerDeriv(y(5),(-params(2)),2);
g2_v(42)=(-(T(1)*T(21)*T(11)*(1-params(2))*exp(y(9))*T(26)));
g2_v(43)=g2_v(42);
g2_v(44)=(-(T(1)*T(21)*(1-params(2))*exp(y(9))*T(10)*T(32)));
g2_v(45)=g2_v(44);
g2_v(46)=(-(T(1)*T(22)));
g2_v(47)=g2_v(46);
g2_v(48)=(-(T(1)*T(12)*T(11)*(1-params(2))*exp(y(9))*T(39)));
g2_v(49)=(-(T(1)*T(12)*(1-params(2))*exp(y(9))*T(26)*T(32)));
g2_v(50)=g2_v(49);
g2_v(51)=(-(T(1)*T(12)*T(11)*(1-params(2))*exp(y(9))*T(26)));
g2_v(52)=g2_v(51);
g2_v(53)=(-(T(1)*T(12)*(1-params(2))*exp(y(9))*T(10)*T(40)));
g2_v(54)=(-(T(1)*T(12)*(1-params(2))*exp(y(9))*T(10)*T(32)));
g2_v(55)=g2_v(54);
g2_v(56)=(-(T(1)*T(13)));
g2_v(57)=(-(T(16)*T(41)));
g2_v(58)=(-(T(23)*T(28)));
g2_v(59)=g2_v(58);
g2_v(60)=(-(T(23)*T(34)));
g2_v(61)=g2_v(60);
g2_v(62)=(-(T(16)*T(23)));
g2_v(63)=g2_v(62);
g2_v(64)=(-(T(17)*T(15)*params(2)*exp(y(9))*getPowerDeriv(y(6),params(2)-1,2)));
g2_v(65)=(-(T(17)*T(27)*T(33)));
g2_v(66)=g2_v(65);
g2_v(67)=(-(T(17)*T(28)));
g2_v(68)=g2_v(67);
g2_v(69)=(-(T(17)*T(14)*getPowerDeriv(y(1),params(2)-1,2)));
g2_v(70)=(-(T(17)*T(34)));
g2_v(71)=g2_v(70);
g2_v(72)=params(12)*exp(y(8));
g2_v(73)=(-(T(16)*T(17)));
g2_v(74)=(-(T(11)*exp(y(9))*T(10)*T(41)));
g2_v(75)=(-(T(23)*T(11)*exp(y(9))*T(26)));
g2_v(76)=g2_v(75);
g2_v(77)=(-(T(23)*exp(y(9))*T(10)*T(32)));
g2_v(78)=g2_v(77);
g2_v(79)=(-(T(11)*exp(y(9))*T(10)*T(23)));
g2_v(80)=g2_v(79);
g2_v(81)=(-(y(1)*(-params(13))+T(17)*T(11)*exp(y(9))*T(39)));
g2_v(82)=(-(T(17)*exp(y(9))*T(26)*T(32)-(params(12)*exp(y(8))+T(4)*2*(y(6)-1))));
g2_v(83)=g2_v(82);
g2_v(84)=(-(y(1)*(-(params(12)*exp(y(8))))));
g2_v(85)=g2_v(84);
g2_v(86)=(-(T(17)*T(11)*exp(y(9))*T(26)));
g2_v(87)=g2_v(86);
g2_v(88)=(-(T(17)*exp(y(9))*T(10)*T(40)));
g2_v(89)=params(12)*exp(y(8))*(y(6)-1);
g2_v(90)=g2_v(89);
g2_v(91)=(-(T(17)*exp(y(9))*T(10)*T(32)));
g2_v(92)=g2_v(91);
g2_v(93)=(-(y(1)*(-(params(12)*exp(y(8))*(y(6)-1)))));
g2_v(94)=(-(T(17)*T(11)*exp(y(9))*T(10)));
g2 = sparse(g2_i,g2_j,g2_v,6,256);
end
