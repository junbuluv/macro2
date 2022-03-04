function [U,V,Uss,Vss,UE,VE]=simulate_U_V(A,B,C,D,N,idx_c,idx_k,idx_n,idx_z,idx_theta,idx_u,N_states,N_controls,N_nodes_total,q_nodes,oo_,sigmaz,sigmat)
    burn_in=500;
    rng('default')
    epsi=[randn(1,N+burn_in)*sigmaz;randn(1,N+burn_in)*sigmat];
    %epsi=zeros(1,N+burn_in);
    %epsi(1)=1.0;
    U=zeros(N_states,N+burn_in);
    V=zeros(N_controls,N+burn_in);
    
    % storage for prime variables need to compute expectations
    UE=zeros(N_states,N+burn_in,N_nodes_total);
    VE=zeros(N_controls,N+burn_in,N_nodes_total);
    
    U(:,1)=B*epsi(:,1);
    V(:,1)=D*epsi(:,1);

    for j=2:N+burn_in
        U(:,j)=A*U(:,j-1)+B*epsi(:,j);
        V(:,j)=C*U(:,j-1)+D*epsi(:,j);
        %%%
        for m=1:N_nodes_total
            UE(:,j,m)=A*U(:,j)+B*q_nodes(m,:)';
            VE(:,j,m)=C*U(:,j)+D*q_nodes(m,:)';
        end
    end

    Uss=[oo_.dr.ys(idx_k); oo_.dr.ys(idx_z); oo_.dr.ys(idx_theta)];
    Vss=[oo_.dr.ys(idx_c); oo_.dr.ys(idx_n); oo_.dr.ys(idx_u)];
    U=U+Uss*ones(1,N+burn_in);
    V=V+Vss*ones(1,N+burn_in);
    
    UE=UE+Uss*ones(1,N+burn_in);
    VE=VE+Vss*ones(1,N+burn_in);
    
     U=U(:,burn_in+1:end);
     V=V(:,burn_in+1:end);
     UE=UE(:,burn_in+1:end,:);
     VE=VE(:,burn_in+1:end,:);
end