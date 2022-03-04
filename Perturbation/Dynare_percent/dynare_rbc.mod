var c n u k z theta;
varexo epsiz epsit;

parameters beta alpha dss B gx gamma rhot rhoz sigmat sigmaz mu phi1 phi2 ;
@#include "params.txt"


model;
    [name='Intertemporal Euler Equation']
    gx*exp(c)^(-gamma) = beta * gx^(1-gamma) * (exp(c(+1))^(-gamma) * ((1- (dss + exp(z(+1))* phi1 * (exp(u(+1))-1) + phi2/2 * (exp(u(+1))-1)^2)) + alpha*exp(theta(+1))*exp(u(+1))^(alpha)*exp(k)^(alpha-1)*exp(n(+1))^(1-alpha)));
    

    [name='Intratemporal Euler Equation']
    B*(1-exp(n))^(-mu)  = exp(c)^(-gamma)*((1-alpha)*exp(theta)*exp(u)^(alpha)*exp(k(-1))^(alpha)*exp(n)^(-alpha));
    

    [name='Capital utilization']
    (exp(z)*phi1 + phi2 * exp((u-1)))  = alpha * exp(theta) * exp(u)^(alpha-1) * exp(k(-1))^(alpha-1) * exp(n)^(1-alpha);
    

    [name='Budget constraint']
    exp(c) + gx*exp(k) = (1- (dss + exp(z)*phi1 * (exp(u) - 1) + phi2/2 * (exp(u)-1)^2))*exp(k(-1)) + exp(theta)*exp(u)^(alpha)*exp(k(-1))^(alpha) * exp(n)^(1-alpha);
    

    [name='Production shock']
    theta = rhot*theta(-1)+epsit;
    

    [name='Utilization shock']
    z = rhoz*z(-1)+epsiz;
end;

initval;
@#include "steady_state.txt"
end;

steady;
check;

shocks;
var epsit; stderr sigmat;
var epsiz; stderr sigmaz;
end;

stoch_simul(order=2, irf= 100, hp_filter = 1600);
