clear all
clc


CPU = cputime;
%% parameter setting
T = 10000; % simulation length
alpha = 1/3; % capital share to output
nss = 1/3; % labor
k2y = 10; % capital to output ratio
i2y = 0.2133; % investment to output ratio
c2y = 1-i2y; % consumption to output ratio
gamma = 1; % consumption risk aversion
mu = 5; % leisure risk aversion
gx = 1.0029; % Labor augmenting rate

sigmaz   = 0.0072; % capital utilization shock standard deviation of disturbance
rhoz = 0.95; % capital utilization shock persistence
sigmat = 0.0096; % production shock standard deviation of disturbance
rhot = 0.95; % production shock persistence

dss = 1-gx + (1- c2y) / k2y; % steady state depreciation rate (u = 1)
beta_s = (gx*k2y)/(c2y-1+alpha+gx*k2y)
phi1 = gx/beta_s - (1-dss); % capital utilization parameter
phi2 = 100; % capital utilization adjustment parameter
B = (1-alpha) * c2y^(-gamma) * k2y^((alpha*(1-gamma))/(1-alpha)) * nss^(-gamma) * (1-nss)^(mu); % leisure utility parameter

%% steady state value calculation
param = struct("beta_s", beta_s, "alpha" ,alpha, "mu", mu, "dss",dss, "gx",gx, "phi1",phi1, "phi2",phi2 , "gamma",gamma,...
    "B", B, "rhot",rhot, "rhoz",rhoz, "sigmat",sigmat,"sigmaz",sigmaz);
init = [0.1,0.1,0.33,1];
options = optimset('Display','iter','MaxFunEvals',1000000,'TolFun',1e-8,'MaxIter',10000);
[ss_val,fval]=fsolve(@(x) steadyPEA(x,param),init,options);


%% steady state
ss_val

css = ss_val(1);
kss = ss_val(2);
nss = ss_val(3);
uss = ss_val(4);



%% Upper bound and lower bound
N(1)=(nss)*0.5;  % the lower bound on n(t)
N(100)=(nss)*2;  % the upper bound on n(t)


%% exogenous shock (production shock, capital utilization shock)
theta = Shocks(T,sigmat,rhot);
z = Shocks(T,sigmaz,rhoz);

%% Initialization 

kp = zeros(T+1,1); % k(t+1)
kp(1) = 0.8*kss;
c = zeros(T,1); % c(t)
n = zeros(T,1); % n(t) 
e = zeros(T-1,1); % Parameterized expectation 
z_p = zeros(T-1,1); % z_p(t)
c_p = zeros(T-1,1); % c_p(t)
u = zeros(T,1); % u(t)
u_p = zeros(T-1,1);% u_p(t)
theta_p = zeros(T-1,1); % theta_p(t) 
n_p = zeros(T-1,1); %n_p(t)
kprime = zeros(T,1);


coef = [log((1-nss)^(-mu)*(gx/beta_s)); 0.001; 0.001; 0.001]; % initial coefficients for labor parametrization

crate    = 0.007;                                   % Speed of moving bounds
criter  	= 1e-6;            				        % Convergence criterion
update  	= .7;            				        % Updating rate 
maxiter  = 100;                                     % Max number of iterations
options = optimset('MaxFunEvals',1000,'TolFun',1e-7,'MaxIter',10000,'Display','off');

%% Quadrature nodes, weight
node_number = 3;
epsi_number = 2;
weight = diag(ones(epsi_number,1));
[n_nodes,epsi_nodes,weight_nodes] = GH_Quadrature(node_number,epsi_number,weight') ;


%% Loop
iteration = 0; % initializing iteration count
dif = Inf;  % norm(zeta - coef)

while (dif > criter)
    up_bound  = kss * 1.5;     % Upper bound
    low_bound = kss * 0.5;         % Lower bound think this process is widening the bounds as iteration goes on. Thus, limiting the interval of kss at the begining.
for t = 1:T
    % Compute the time series
    % Parameterization of n(t)
    n(t) = 1- (beta_s/gx * exp(coef(1) + coef(2) * log(kp(t)) +  coef(3) * log(theta(t)) + coef(4) * log(z(t))))^(-1/mu);
    n(t)=n(t)*(n(t)>=N(1))*(n(t)<=N(100))+N(1)*(n(t)<N(1))+N(100)*(n(t)>N(100)); %making n(t) not go over the bounds 
    % calculate u(t) using FOC
    u(t)=((alpha*theta(t)*kp(t)^(alpha-1)*n(t)^(1-alpha) - z(t)*phi1)/phi2) +1;
    %calculate c(t) using FOC
    c(t) = (B*(1-n(t))^(-mu) / ((1-alpha)*theta(t)*u(t)*kp(t)^(alpha)*n(t)^(-alpha)))^(-1/gamma);
    %calculate kp using Budget constraint
    kp(t+1) = (1/gx)*((1- (dss + phi1 * (u(t)-1) + phi2/2 * (u(t)-1)^2))*kp(t) + theta(t)*u(t)*kp(t)^(alpha)*n(t)^(1-alpha) - c(t));

    %calculate u(t) counterpart
    if kp(t+1) > up_bound  
          kp(t+1) = up_bound; 
   elseif kp(t+1) < low_bound
          kp(t+1) = low_bound;
    end

end



%% Expectation part calculation
% z_prime and theta_prime
z_p = z.^(rhoz) * exp(sigmaz * epsi_nodes(:,2))';
theta_p = theta.^(rhot) * exp(sigmat * epsi_nodes(:,1))';
k = repmat(kp(1:end-1),1,n_nodes);
kprime = repmat(kp(2:end),1,n_nodes);

%labor_prime - parametrization
n_p = (1- (beta_s/gx .* exp(coef(1) + coef(2).* log(kprime) +  coef(3) .* log(theta_p) + coef(4).* log(z_p))).^(-1/mu)) * weight_nodes;
%capital utilization prime - FOC
u_p = (((alpha*theta_p.*kprime.^(alpha-1).*n_p.^(1-alpha) - z_p.*phi1)./phi2) + 1 ) * weight_nodes;
%consumption prime - FOC
c_p = ((B.*(1-n_p).^(-mu) ./ ((1-alpha)*theta_p.*u_p.*kprime.^(alpha).*n_p.^(-alpha))).^(-1/gamma)) * weight_nodes;
% Gauss-Hermite Quadrature
z_p = z_p * weight_nodes;
theta_p = theta_p * weight_nodes;
% recover capital and capital prime Quadrature
k = k*weight_nodes;
kprime = kprime * weight_nodes;

e = (1-n_p).^(-mu) ./ (theta.* u.* k.^(alpha).* n.^(-alpha)) .* (theta_p.*u_p.*kprime.^(alpha).*n_p.^(-alpha)) .*...
    ((1- (dss + z_p.* phi1 .* (u_p - 1) + phi2/2 .* (u_p - 1).^2)) + alpha .* theta_p .* u_p .* kprime.^(alpha-1) .* n_p.^(1-alpha));

X = [ones(T,1), log(kp(1:end-1)), log(theta(1:end)), log(z(1:end))];
zeta = nlinfit(X,e,'object',coef);


dif = norm(coef - zeta)
if dif > criter
    coef = update*zeta + (1-update)*coef;
else
    fprintf("Optimal Coefficient for parametrized function is %.3f \n", coef);
end
iteration = iteration+1;

end



% 7.0 Computational time
% ----------------------
CPU = cputime-CPU


time=(1:1:T);                         
title('Time series solution');
subplot(4,1,1);
plot (time,kp(1:T,1)), xlabel('t'), ylabel('Capital')
subplot(4,1,2);
plot (time,c), xlabel('t'), ylabel('Consumption')
subplot(4,1,3);
plot (time,n), xlabel('t'), ylabel('labor')
subplot(4,1,4);
plot (time,u), xlabel('t'), ylabel('capital utilization')

