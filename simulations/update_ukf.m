function [xe,Pe,lm_seq] = update_ukf(xe,Pe,lm_seq,z,R)


lenz = length(find(z(3,:)>0));

% % update
Idx = [];
z_mes = [];
for i= 1:lenz
    % data association (based on landmark id)
    % TODO: use nearest neighbor to do the job
    is_exist = ~(lm_seq - z(3,i));
    idx = find(is_exist);
    if ~isempty(idx)
        Idx = [Idx idx];
        z_mes = [z_mes;z(1:2,i)];
    end
end

% update: already in the state vecor
if ~isempty(Idx)
    RR = kron(eye(length(Idx)),R(1:2,1:2));
    L = length(Pe);
    [c,Ws,Wc] = unscented_transform_parameters(L,0.2);
    P_root = chol(Pe)'; % make sure P is positive definite
    dX = [ zeros(L,1), c*P_root(:,1:L), -c*P_root(:,1:L) ];
    
    Z = zeros(2*length(Idx), 2*L+1 );
    for k = 1:2*L+1
        dX_k = dX(:,k);
        x = xe+dX_k;
        for j = 1:length(Idx)
            Z(2*j-1:2*j,k) = measurement_model_ukf(x,Idx(j));
        end
    end
    z_hat = Ws(1)*Z(:,1)+Ws(2)*sum(Z(:,2:end),2);
    
    Pxz = (dX)*diag(Wc)*(Z-z_hat*ones( 1, 2*L+1 ))';
    H = Pxz' / Pe;
    
    Pzz = (Z-z_hat*ones( 1, 2*L+1 ))*diag(Wc)*(Z-z_hat*ones( 1, 2*L+1 ))';
    Pzz = Pzz + RR;
    Kgain = Pe*H'/Pzz;
    
    dx = Kgain*(z_mes-z_hat);
    Pe = Pe-Kgain*Pzz*Kgain';
    
    xe = xe + dx;
    
    
end




% % augment
for i= 1:lenz
    % data association (known)
    is_exist = ~(lm_seq - z(3,i));
    Idx = find(is_exist);
    
    lenx= size(xe,1);
    ii = 2*i + (-1:0);
    
    % add the new landmark into the state vector
    if isempty(Idx)
        lm_seq = [lm_seq; z(3,i)];
        
        % augment state
        k_xL = z(1:2,i);
        
        C = [cos(xe(3)) -sin(xe(3));  sin(xe(3)) cos(xe(3)) ];
        x_L = xe(1:2,1) + C*k_xL;
        xe = [xe; x_L];
        
        % jacobians
        J = [0 -1; 1 0];
        C = [cos(xe(3))  -sin(xe(3));   sin(xe(3))   cos(xe(3)) ];
        
        HR = - C'*[eye(2)  J*(x_L-xe(1:2,1))];
        HL = C';
        
        % augment covariance
        rng= lenx+1:lenx+2;
        Pe(rng,rng)= inv(HL)*HR*Pe(1:3,1:3)*HR'*inv(HL)' + inv(HL)*R(ii,ii)*inv(HL)'; % landmark cov
        Pe(rng,1:3)= -inv(HL)*HR*Pe(1:3,1:3); % landmark-robot xcorr
        Pe(1:3,rng)= Pe(rng,1:3)';
        if lenx>3
            rnm= 4:lenx;
            Pe(rng,rnm)= -inv(HL)*HR*Pe(1:3,rnm);
            Pe(rnm,rng)= Pe(rng,rnm)';
        end
        Pe = 0.5*(Pe+Pe');
    end
end

