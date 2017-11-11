figure;
plot(xRest_std(1,:,1),xRest_std(2,:,1)','b-','Linewidth',2);
hold on;
plot(xRest_ocekf(1,:,1),xRest_ocekf(2,:,1)','k-','Linewidth',2);
plot(xRest_iekf(1,:,1),xRest_iekf(2,:,1)','c-','Linewidth',2);
plot(xRest_rc(1,:,1),xRest_rc(2,:,1)','g-','Linewidth',2);
plot(xR_true(1,:),xR_true(2,:)','y-','Linewidth',2);
% comparaison_gps
% plot(x_GPS(1,:),x_GPS(2,:)','r-','Linewidth',1);
