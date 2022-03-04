





%% get policy function from Dynare


idx_c = 1;
idx_n = 2;
idx_u = 3;
idx_k = 4;
idx_z = 5;
idx_theta = 6;

[A,B,C,D] = get_A_B_C_D(idx_c,idx_k,idx_n,idx_u,idx_z,idx_theta,oo_);

% move into a function
% number of nodes we want for each innovation
N_nodes_individual=3;
% number of innovations
N_innov=2;
% variance covariance matrix
vcv_mat=1.0;
[N_nodes_total,qnodes,qweights] = GH_Quadrature(N_nodes_individual,N_innov,vcv_mat);


N=1000;
N_states=3;
N_controls=3;
[U,V,Uss,Vss,UE,VE]= simulate_U_V(A,B,C,D,N,idx_c,idx_k,idx_n,idx_z,idx_theta,idx_u,N_states,N_controls,N_nodes_total,qnodes,oo_);


idx_c_v=1;
idx_n_v=2;
idx_u_v=3;
idx_k_u=1;
idx_z_u=2;
idx_theta_u=3;


% Construct residuals
sim_len=length(U(1,:));
E=zeros(sim_len,N_nodes_total);

for j=2:sim_len
    for m=1:N_nodes_total
       E(j,m)=((param.beta*VE(idx_c_v,j,m).^(-param.gamma).*(param.alpha*exp(UE(idx_z_u,j,m)).*(U(idx_k_u,j)).^(param.alpha-1.0)+(1.0-param.dss))));
    end
end

% apply weights
E=E*qweights;

% compute errors
EERES=E(2:end).^(-1.0/param.gamma)./(V(idx_c_v,2:end)')-1.0;

max_ee_res=log10(max(abs(EERES)));

fprintf("Max Euler equation residual is %.3f \n",max_ee_res);


