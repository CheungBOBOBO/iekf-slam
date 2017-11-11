function [zhat] = measurement_model_ukf(xe,idx)
%
% relative position meausrement model
%
fpos = 3+idx*2-1;
C = [cos(xe(3,1)) -sin(xe(3,1));
    sin(xe(3,1)) cos(xe(3,1))];

k_xL = C'*(xe(fpos:fpos+1)-xe(1:2));
zhat = k_xL;
end

