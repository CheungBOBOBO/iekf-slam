close all

% % setup 
for nr= 1%1:nRuns
    figure
    plot(xR_true(1,:,nr), xR_true(2,:,nr), 'b-.','Linewidth',1), hold on
    plot(xL_true_fixed(1,:), xL_true_fixed(2,:), 'b*','Linewidth',1)
    title('Robot trajectory and landmarks')
    legend('True Trajectory','True Landmarks')
    xlabel('x (m)','FontWeight','bold'), ylabel('y (m)','FontWeight','bold')
    axis equal
end

% % robot % %

start= 1;  incr= 5;

% Robot NEES
figure, hold on
plot([start:incr:nSteps],neesR_avg_id(start:incr:end)','g:','Linewidth',2);
plot([start:incr:nSteps],neesR_avg_std(start:incr:end)','b-','Linewidth',2);
plot([start:incr:nSteps],neesR_avg_iekf(start:incr:end)','m-.','Linewidth',2);
plot([start:incr:nSteps],neesR_avg_ocekf(start:incr:end)','k-','Linewidth',2);
plot([start:incr:nSteps],neesR_avg_rc(start:incr:end)','c-','Linewidth',2);
plot([start:incr:nSteps],neesR_avg_ukf(start:incr:end)','y-','Linewidth',2);
xlabel('Time (sec)','FontWeight','bold'), ylabel('Robot pose NEES','FontWeight','bold')
legend( 'Ideal-EKF','Std-EKF','IEKF','OC-EKF','Robocentric','Std-UKF')


% Robot RMSE
figure
subplot(2,1,1), hold on
plot([start:incr:nSteps],rmsRp_avg_id(start:incr:end)','g:','Linewidth',2);
plot([start:incr:nSteps],rmsRp_avg_std(start:incr:end)','b-','Linewidth',2);
plot([start:incr:nSteps],rmsRp_avg_iekf(start:incr:end)','m-.','Linewidth',2);
plot([start:incr:nSteps],rmsRp_avg_ocekf(start:incr:end)','k-','Linewidth',2);
plot([start:incr:nSteps],rmsRp_avg_rc(start:incr:end)','c-','Linewidth',2);
plot([start:incr:nSteps],rmsRp_avg_ukf(start:incr:end)','y-','Linewidth',2);
ylabel('Position RMSE (m)','FontWeight','bold')
legend( 'Ideal-EKF','Std-EKF','IEKF','OC-EKF','Robocentric','Std-UKF')
subplot(2,1,2), hold on
plot([start:incr:nSteps],rmsRth_avg_id(start:incr:end)','g:','Linewidth',2);
plot([start:incr:nSteps],rmsRth_avg_std(start:incr:end)','b-','Linewidth',2);
plot([start:incr:nSteps],rmsRth_avg_iekf(start:incr:end)','m-.','Linewidth',2);
plot([start:incr:nSteps],rmsRth_avg_ocekf(start:incr:end)','k-','Linewidth',2);
plot([start:incr:nSteps],rmsRth_avg_rc(start:incr:end)','c-','Linewidth',2);
plot([start:incr:nSteps],rmsRth_avg_ukf(start:incr:end)','y-','Linewidth',2);
xlabel('Time (sec)','FontWeight','bold'),
ylabel('Heading RMSE (rad)','FontWeight','bold')


NEES_Robot = [mean(neesR_avg_id(5:end)),mean(neesR_avg_std(5:end)), mean(neesR_avg_iekf(5:end)), mean(neesR_avg_ocekf(5:end)) ...
    , mean(neesR_avg_ukf(5:end)), mean(neesR_avg_rc(5:end))];
RMS_Position = [mean(rmsRp_avg_id(5:end)),mean(rmsRp_avg_std(5:end)), mean(rmsRp_avg_iekf(5:end)), ...
    mean(rmsRp_avg_ocekf(5:end)), mean(rmsRp_avg_ukf(5:end)), mean(rmsRp_avg_rc(5:end)) ];
RMS_Heading = [mean(rmsRth_avg_id(5:end)),mean(rmsRth_avg_std(5:end)), mean(rmsRth_avg_iekf(5:end)), ...
    mean(rmsRth_avg_ocekf(5:end)), mean(rmsRth_avg_ukf(5:end)), mean(rmsRth_avg_rc(5:end))];

disp(NEES_Robot);
disp(RMS_Position);
disp(RMS_Heading);

