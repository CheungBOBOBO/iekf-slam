function [c,Ws,Wc] = unscented_transform_parameters(L,alpha)
kappa  = 3-L;
beta   = 2;
lambda = alpha^2*(L+kappa)-L;
Ws     = [ lambda/(L+lambda), 1/(2*(L+lambda))*ones(1,2*L) ];
Wc     = [ lambda/(L+lambda)+(1-alpha^2+beta), 1/(2*(L+lambda))*ones(1,2*L) ];
c      = (L+lambda)^0.5;
end

