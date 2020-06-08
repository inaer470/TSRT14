function [xhat, meas] = filterTemplate(fname, calAcc, calGyr, calMag)
% FILTERTEMPLATE  Filter template
%
% This is a template function for how to collect and filter data
% sent from a smartphone live.  Calibration data for the
% accelerometer, gyroscope and magnetometer assumed available as
% structs with fields m (mean) and R (variance).  If fname is provided and
% not empty, the log data in fname is used instead of streamed data (the
% timing is not accurate).
%
% The function returns xhat as an array of structs comprising t
% (timestamp), x (state), and P (state covariance) for each
% timestamp, and meas an array of structs comprising t (timestamp),
% acc (accelerometer measurements), gyr (gyroscope measurements),
% mag (magnetometer measurements), and orint (orientation quaternions
% from the phone).  Measurements not availabe are marked with NaNs.
%
% As you implement your own orientation estimate, it will be
% visualized in a simple illustration.  If the orientation estimate
% is checked in the Sensor Fusion app, it will be displayed in a
% separate view.

  %% Setup necessary infrastructure
  import('se.hendeby.sensordata.*');  % Used to receive data.

  %% Filter settings
  t0 = [];  % Initial time (initialize on first data received)
  nx = 4;
  % Add your filter settings here.

  % Current filter state.
  x = [1; 0; 0 ;0];
  P = eye(nx, nx);

  % Saved filter states.
  xhat = struct('t', zeros(1, 0),...
                'x', zeros(nx, 0),...
                'P', zeros(nx, nx, 0));

  meas = struct('t', zeros(1, 0),...
                'acc', zeros(3, 0),...
                'gyr', zeros(3, 0),...
                'mag', zeros(3, 0),...
                'orient', zeros(4, 0));
  try
    %% Create data link
    if nargin == 0 || isempty(fname)
      server = StreamSensorDataReader(3400);
    else
      server = FileSensorDataReader(fname);
    end
    % Makes sure to resources are returned.
    sentinel = onCleanup(@() server.stop());

    server.start();  % Start data reception.
  catch e
    fprintf(['Unsuccessful connecting to client!\n' ...
      'Make sure to start streaming from the phone *after* '...
             'running this function!']);
    return;
  end

  % Used for visualization.
  figure(1);
  subplot(1, 2, 1);
  ownView = OrientationView('Own filter', gca);  % Used for visualization.
  ownView.activateKeyboardCallback;
  googleView = [];
  counter = 0;  % Used to throttle the displayed frame rate.

  %% Filter loop
  % Repeat while data is available and q hasn't been pressed
  while server.status() && ~ownView.quit
    % Get the next measurement set, assume all measurements
    % within the next 5 ms are concurrent (suitable for sampling
    % in 100Hz).
    data = server.getNext(5);

    if isnan(data(1))  % No new data received
      continue;
    end
    t = data(1)/1000;  % Extract current time

    if isempty(t0)  % Initialize t0
      t0 = t;
    end

    gyr = data(1, 5:7)';
    if ~any(isnan(gyr))  % Gyro measurements are available.
      % Do something
      
    end

    acc = data(1, 2:4)';
    if ~any(isnan(acc))  % Acc measurements are available.
      % Do something
    end

    mag = data(1, 8:10)';
    if ~any(isnan(mag))  % Mag measurements are available.
      % Do something
    end

    orientation = data(1, 18:21)';  % Google's orientation estimate.

    % Visualize result
    if rem(counter, 10) == 0
      setOrientation(ownView, x(1:4));
      title(ownView, 'OWN', 'FontSize', 16);
      if ~any(isnan(orientation))
        if isempty(googleView)
          subplot(1, 2, 2);
          % Used for visualization.
          googleView = OrientationView('Google filter', gca);
        end
        setOrientation(googleView, orientation);
        title(googleView, 'GOOGLE', 'FontSize', 16);
      end
    end
    counter = counter + 1;

    % Save estimates
    xhat.x(:, end+1) = x;
    xhat.P(:, :, end+1) = P;
    xhat.t(end+1) = t - t0;

    meas.t(end+1) = t - t0;
    meas.acc(:, end+1) = acc;
    meas.gyr(:, end+1) = gyr;
    meas.mag(:, end+1) = mag;
    meas.orient(:, end+1) = orientation;
  end
%    mean and histogram
  mean_acc = mean(meas.acc(:, ~any(isnan(meas.acc), 1)), 2)
  mean_gyr = mean(meas.gyr(:, ~any(isnan(meas.gyr), 1)), 2)
  mean_mag = mean(meas.mag(:, ~any(isnan(meas.mag), 1)), 2)
  
  var_acc1 = var(meas.acc(1, ~any(isnan(meas.acc), 1)),1);
  var_acc2 = var(meas.acc(2, ~any(isnan(meas.acc), 1)),1);
  var_acc3 = var(meas.acc(3, ~any(isnan(meas.acc), 1)),1);
  
  var_acc = diag([var_acc1, var_acc2, var_acc3])
  
  var_gyr1 = var(meas.gyr(1, ~any(isnan(meas.gyr), 1)),1);
  var_gyr2 = var(meas.gyr(2, ~any(isnan(meas.gyr), 1)),1);
  var_gyr3 = var(meas.gyr(3, ~any(isnan(meas.gyr), 1)),1);
  
  var_gyr = diag([var_gyr1, var_gyr2, var_gyr3])
  
  var_mag1 = var(meas.mag(1, ~any(isnan(meas.mag), 1)),1);
  var_mag2 = var(meas.mag(2, ~any(isnan(meas.mag), 1)),1);
  var_mag3 = var(meas.mag(3, ~any(isnan(meas.mag), 1)),1);
  
  var_mag = diag([var_mag1, var_mag2, var_mag3])
  
  
  figure(2)
  subplot(3,3,1)
  histfit(meas.acc(1, ~any(isnan(meas.acc), 1)));
  title('acceletation x')
  subplot(3,3,2)
  histfit(meas.acc(2, ~any(isnan(meas.acc), 1)));
  title('acceletation y')
  subplot(3,3,3)
  histfit(meas.acc(3, ~any(isnan(meas.acc), 1)));
  title('acceletation z')
  
  subplot(3,3,4)
  histfit(meas.gyr(1, ~any(isnan(meas.gyr), 1)));
  title('gyroscope x')
  subplot(3,3,5)
  histfit(meas.gyr(2, ~any(isnan(meas.gyr), 1)));
  title('gyroscope y')
  subplot(3,3,6)
  histfit(meas.gyr(3, ~any(isnan(meas.gyr), 1)));
  title('gyroscope z')
  
  subplot(3,3,7)
  histfit(meas.mag(1, ~any(isnan(meas.mag), 1)));
  title('magnetic field x')
  subplot(3,3,8)
  histfit(meas.mag(2, ~any(isnan(meas.mag), 1)));
  title('magnetic field y')
  subplot(3,3,9)
  histfit(meas.mag(3, ~any(isnan(meas.mag), 1)));
  title('magnetic field z')
% figure(2)
% subplot(3,1,1)
% plot(meas.acc(1, ~any(isnan(meas.acc), 1)))
% hold on
% plot(meas.acc(2, ~any(isnan(meas.acc), 1)))
% plot(meas.acc(3, ~any(isnan(meas.acc), 1)))
% title('accelerometer')
% legend('x','y','z')
% 
% subplot(3,1,2)
% plot(meas.gyr(1, ~any(isnan(meas.gyr), 1)))
% hold on
% plot(meas.gyr(2, ~any(isnan(meas.gyr), 1)))
% plot(meas.gyr(3, ~any(isnan(meas.gyr), 1)))
% title('gyroscope')
% legend('x','y','z')
% 
% subplot(3,1,3)
% plot(meas.mag(1, ~any(isnan(meas.mag), 1)))
% hold on
% plot(meas.mag(2, ~any(isnan(meas.mag), 1)))
% plot(meas.mag(3, ~any(isnan(meas.mag), 1)))
% title('magnetic field')
% legend('x','y','z')
end
