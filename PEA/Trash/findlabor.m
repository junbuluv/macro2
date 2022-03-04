function val = findlabor(n,mp,c,k,u,theta)

val = mp.B*(1-n)^(-mp.mu) - c^(-mp.gamma) *((1-mp.alpha) * theta * (u*k)^(mp.alpha)*n^(-mp.alpha)) ;
end