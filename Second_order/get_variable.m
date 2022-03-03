function [c_p,n_p,u_p,z_p,theta_p] = get_variable(sim_len,N_nodes_total,VE,UE,idx_c_v,idx_n_v,idx_u_v,idx_z_u,idx_theta_u, qweights)

c_p = zeros(sim_len,N_nodes_total);
n_p = zeros(sim_len,N_nodes_total);
u_p = zeros(sim_len,N_nodes_total);
z_p = zeros(sim_len,N_nodes_total);
theta_p = zeros(sim_len,N_nodes_total);

for j = 2:sim_len
    for m = 1:N_nodes_total
        c_p(j,m) = VE(idx_c_v,j,m);
        n_p(j,m) = VE(idx_n_v,j,m);
        u_p(j,m) = VE(idx_u_v,j,m);
        z_p(j,m) = UE(idx_z_u,j,m);
        theta_p(j,m) = UE(idx_theta_u,j,m);
    end
end


c_p = c_p * qweights;
n_p = n_p * qweights;
u_p = u_p * qweights;
z_p = z_p * qweights;
theta_p = theta_p * qweights;


end