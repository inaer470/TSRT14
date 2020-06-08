clear all
%% Task 1
load('/Users/Ina/Desktop/TSRT14 -LAB/Lab1/SFlab1-data/calibration.mat');
% Remove the second mic
nr_of_mics = 7;
tphat(:,2) = [];

% Scale from seconds to meter
rhat = tphat*343; 


%%

for i = 1:nr_of_mics
    % Calculate measurement error, bias, variance and 
    % standard deviation for each microphone
    e(:,i) = rhat(:,i) - mean(rhat, 2);
    bias(i) = mean(e(:,i));
    variance(i) = var(e(:,i));
    standDev(i) = std(e(:,i));
    
    % Plot the histogram
%     subplot(2, 4, i);
%     histfit(e(:,i))
%     title(['Microphone ', num2str(i)]);

end
% sgtitle('Histogram of measurement error and normal distribution');

%% task 3
clear tphat;
clear rhat;

load('/Users/Ina/Desktop/TSRT14 -LAB/Lab1/SFlab1-data/setup1.mat');
% Remove the second mic
tphat(:,2) = [];

% Scale from seconds to meter
rhat = tphat*343;
%Scale from centimeters to meter
mic_locations = mic_locations'/100;
mic_locations(:,2) = [];
 x0 = x0/100;

% New locations used in Sensitivity analysis
%  mic_locations =  [0       0.3900    0.8100    1.2220    1.2220    0.8100    0.39000;
%                  0.6100          0         0    0.2900    0.6100    0.9910    0.9910];


% Create the first sensormod object
sm_TDOA1 = exsensor('tdoa1', nr_of_mics);
sm_TDOA1.th = mic_locations(:);
sm_TDOA1.x0 = [x0, 0];
sm_TDOA1.pe = diag(variance);

% Create the second sensormod object
h = inline('[sqrt((x(1,:)-th(13)).^2+(x(2,:)-th(14)).^2) - sqrt((x(1,:)-th(1)).^2+(x(2,:)-th(2)).^2);sqrt((x(1,:)-th(13)).^2+(x(2,:)-th(14)).^2) - sqrt((x(1,:)-th(3)).^2+(x(2,:)-th(4)).^2);sqrt((x(1,:)-th(13)).^2+(x(2,:)-th(14)).^2) - sqrt((x(1,:)-th(5)).^2+(x(2,:)-th(6)).^2);sqrt((x(1,:)-th(13)).^2+(x(2,:)-th(14)).^2) - sqrt((x(1,:)-th(7)).^2+(x(2,:)-th(8)).^2);sqrt((x(1,:)-th(13)).^2+(x(2,:)-th(14)).^2) - sqrt((x(1,:)-th(9)).^2+(x(2,:)-th(10)).^2);sqrt((x(1,:)-th(13)).^2+(x(2,:)-th(14)).^2) - sqrt((x(1,:)-th(11)).^2+(x(2,:)-th(12)).^2)]', 't', 'x', 'u', 'th');  
sm_TDOA2 = sensormod(h, [2 0 6 14]);
sm_TDOA2.th = mic_locations(:);
sm_TDOA2.x0 = x0';
sm_TDOA2.pe = get_R(variance);

%%

%   figure(1)
%   plot(sm_TDOA1)
%    scatter(mic_locations(1,:), mic_locations(2,:),70, '*');
%   hold on
%     scatter(mic_locations2(1,:), mic_locations2(2,:),70, '*r');
%   xlabel('Meter (m)')
%   ylabel('Meter (m)')
%   title('Cramér-Rao Lower Bound for TDOA1 with setup 2')
%   hold on
%   figure(2)
%   plot(sm_TDOA2)
%   title('Cramér-Rao Lower Bound for TDOA2 with setup 2')
%   hold on
%   xlabel('Meter (m)')
%   ylabel('Meter (m)')

%% Task 4 crlb
% Plot the CRLB for TDOA1 grid and display the RMSE
%  figure(1)
% crlb2(sm_TDOA1, [], 0:.1:1.3, 0:.1:1, [1,2], 'rmse');
% axis([0 1 0 1])
%   xlabel('Meter (m)')
%   ylabel('Meter (m)')
% % To plot 3D graph
% % [cx, X1, X2] = crlb2(sm_TDOA1, [], 0:.1:1.3, 0:.1:1, [1,2], 'rmse');
% % surf(X1, X2, cx);
% 
% % Plot the CRLB for TDOA1 grid and display the RMSE
% figure(2)
% crlb2(sm_TDOA2, [], 0:.1:1.3, 0:.1:1, [1,2], 'rmse');
% 
% axis([0 1 0 1])
%   xlabel('Meter (m)')
%   ylabel('Meter (m)')

%% task 5 b) TDOA1
% 
% % Get measurements for TDOA1 using function 
y = y_TDOA1(rhat);

% Set the start position to a position in the middle
sm_TDOA1.x0 = [x0, min(rhat(1,:))-0.002];

for i = 1:length(rhat)-1
    % , 'maxiter', 1000, 'ctol', 1e-7
    % Use estimate to estimate the position
    [shat, xhat] = estimate(sm_TDOA1, y(i,:));
    
    % Store the estimated position
    estimated_pos(:,i) = shat.x0(1:2);

    % Update r0
    sm_TDOA1.x0 = [shat.x0(1:2); min(rhat(i+1,:))-0.002];
%    
%     
%     plot(shat, 'conf', 90)
%     hold on
end
% 
% SFlabCompEstimGroundTruth(estimated_pos, mic_locations);
% set(gca, 'xTickLabel', []);
% set(gca, 'yTickLabel', []);

%% task 5 d) TDOA2

% % Get measurements for TDOA2 using function 
% y = y_TDOA2(rhat, mic_locations);
% 
% % Set the start position to a position in the middle
% sm_TDOA2.x0 = [x0];
% 
% for i = 1:length(rhat)
%     % Use estimate to estimate the position
%     [xhat2, shat] = estimate(sm_TDOA2, y(i,:));
%     
%     %Store the estimated position
%     estimated_pos(:,i) = xhat2.x0(1:2);
%     sm_TDOA2.x0 = [xhat2.x0(1:2)];
%     %plot(xhat2, 'conf', 90)
%     %hold on
% end
% 
% SFlabCompEstimGroundTruth(estimated_pos, mic_locations);
% set(gca, 'xTickLabel', []);
% set(gca, 'yTickLabel', []);

%% Task 6 a 
% % Load workspace from task 5 b
% load('artificial_measurments_sensitive.mat') 
% 
% z = estimated_pos;
% 
% % Change motion model (cv2D or ca2D)
% m = exlti('cv2D');
% %m = exlti('ca2D');
% 
% % Improving the filter
% m.Q = 1 * m.Q;
% 
% y = sig(z');
% % Apply the Kalman filter
% xhat_TDOA1 = kalman(m, y);
% 
% %  xplot2(xhat_TDOA1, z, 'conf', 99);
% % Plot the confidence bounds
% % xplot(xhat_TDOA1, sig(z).y, 'conf', 99);
% % 
% % Plot the trajectory 
% SFlabCompEstimGroundTruth(xhat_TDOA1.y', mic_locations);
% title('Estimated trajectory for CV with KF, disturbed positions')
% set(gca, 'xTickLabel', []);
% set(gca, 'yTickLabel', []);

%% Task 6 b
% 
% 
% y = y_TDOA2(rhat, mic_locations);
% 
% % Change motion model (cv2D or ca2D)
% m = exmotion('cv2D');
% %m = exmotion('ca2D');
% 
% 
% % Combine the motion and sensor models
% m1 = addsensor(m, sm_TDOA2);
% 
% % Improving the filter
% m1.pv = 0.1 * m1.pv;
% m1.pe = 10* m1.pe;
% %m1.x0 = [x0 0.1 0.1 0 0]; for ca
% m1.x0 = [x0 0.1 0.1];
% 
% % Apply the extended Kalman filter
% xhat_TDOA2 = ekf(m1, y);
% 
% 
% % xplot(xhat_TDOA1, sig(z).y, 'conf', 99);
% 
% 
% % Plot the trajectory 
% SFlabCompEstimGroundTruth(xhat_TDOA2.x(:,1:2)', mic_locations);
% set(gca, 'xTickLabel', []);
% set(gca, 'yTickLabel', []);
% 
% % Plot trajextories for both TDOA1 and TDOA2 
% xplot2(xhat_TDOA2(1:50,1:2), xhat_TDOA1(1:50,:));
% 
% title('Estimated trajectory using motion model CV')
% xlabel('Meter (m)')
% ylabel('Meter (m)')
