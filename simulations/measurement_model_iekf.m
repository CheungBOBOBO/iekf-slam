function [zhat, H] = measurement_model_iekf(xe,idx)
%
% relative position meausrement model
%



fpos = 3+idx*2-1;


C = [cos(xe(3,1)) -sin(xe(3,1));
    sin(xe(3,1)) cos(xe(3,1))];



zhat = rot(xe(3))'*(xe(fpos:fpos+1)-xe(1:2));


H = zeros(2,size(xe,1));
H(:,1:2) = -rot(xe(3))';
H(:,fpos:fpos+1) =  rot(xe(3))';

end
