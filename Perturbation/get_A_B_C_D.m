function [A,B,C,D]=get_A_B_C_D_E(idx_c,idx_k,idx_n,idx_u,idx_z,idx_theta,oo_)

    A=[oo_.dr.ghx(oo_.dr.inv_order_var(idx_k),:);...
       oo_.dr.ghx(oo_.dr.inv_order_var(idx_z),:);
       oo_.dr.ghx(oo_.dr.inv_order_var(idx_theta),:)];

    B=[oo_.dr.ghu(oo_.dr.inv_order_var(idx_k),:);...
       oo_.dr.ghu(oo_.dr.inv_order_var(idx_z),:);
       oo_.dr.ghu(oo_.dr.inv_order_var(idx_theta),:)];

    C=[oo_.dr.ghx(oo_.dr.inv_order_var(idx_c),:);...
       oo_.dr.ghx(oo_.dr.inv_order_var(idx_n),:);
       oo_.dr.ghx(oo_.dr.inv_order_var(idx_u),:)];

    D=[oo_.dr.ghu(oo_.dr.inv_order_var(idx_c),:);...
       oo_.dr.ghu(oo_.dr.inv_order_var(idx_n),:);
       oo_.dr.ghu(oo_.dr.inv_order_var(idx_u),:)];
end
