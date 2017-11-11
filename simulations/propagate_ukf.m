function [xe, Pe] = propagate_ukf(xe,Pe,dt,v_m,omega_m,sigma_v,sigma_w)
Q = [sigma_v^2  0;   0    sigma_w^2];
% % propagate state
xe(1:3,1) = [ xe(1) + v_m*dt*cos(xe(3,1));
    xe(2) + v_m*dt*sin(xe(3,1));
    pi_to_pi(xe(3) + omega_m*dt) ];


idxP = 1:3;
Pxx = Pe(idxP,idxP);
P_root = chol(blkdiag(Pxx,Q))'; % make sure P is positive definite
L = length(P_root);
[c,Ws,Wc] = unscented_transform_parameters(L,0.2);
dX = [ zeros(L,1), c*P_root(:,1:L), -c*P_root(:,1:L) ];

Z = zeros(3, 2*L+1 );
for k = 1:2*L+1
    dX_k = dX(:,k);
    Z(:,k) =  [ dX_k(1) + (v_m+dX_k(4))*dt*cos(xe(3,1)+dX_k(3));
    dX_k(2) + (v_m+dX_k(4))*dt*sin(xe(3,1)+dX_k(3));
    dX_k(3) + (omega_m+dX_k(5))*dt ];
end
z_hat = Ws(1)*Z(:,1)+Ws(2)*sum(Z(:,2:end),2);
Pxz = (dX)*diag(Wc)*(Z-z_hat*ones( 1, 2*L+1 ))';
A = Pxz(1:3,1:3)' / Pxx;
B = zeros(3,length(Pe));
B(:,1:3) = A;
Pzz = (Z-z_hat*ones( 1, 2*L+1 ))*diag(Wc)*(Z-z_hat*ones( 1, 2*L+1 ))';

Pe = [Pzz B*Pe(:,4:end);Pe(4:end,:)*B' Pe(4:end,4:end)];
end
