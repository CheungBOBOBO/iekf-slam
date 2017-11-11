function xL = gen_map(nL,v,omega,rmin,rmax, nSteps,dt)

% radius
r = v/omega;

for i=1:nL
    rho = r + rmin + 2;
    th =   4/4*2*pi*i/nL;
    [x,y] = pol2cart(th,rho);
    xL(:,i) = [x;y+r]; %shift y w/ r since robot starts at (0,0)
end

end

