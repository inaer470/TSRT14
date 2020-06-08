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
  
% mean_acc =[0.0187; 0.0089; 9.8318];
% 
% mean_gyr =[-0.0006; 0.0009; -0.0011];
% 
% mean_mag =[ -7.5403; 1.0386; -52.0281];
mean_acc =[-0.0333;0.0178;9.8410];
mean_gyr =[0.0000;0.0007;-0.0014];
mean_mag =[ -9.3386; 8.0742; -51.9027];

var_acc =1.0e-03 * [0.0547         0         0;
                     0    0.0406         0;
                        0         0    0.1655];


var_gyr = 1.0e-04 *[0.0619         0         0;
                        0    0.2730        0;
                         0         0    0.0027];


var_mag = [0.1267         0         0;
             0    0.1530         0;
             0         0    0.1509]

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
         [x,P] = tu_qw(x, P, gyr-mean_gyr, meas.t(end)-meas.t(end-1), var_gyr);
    end

    acc = data(1, 2:4)';
    if ~any(isnan(acc))  % Acc measurements are available.
      % Do something
           if (norm(acc,2) < 11 && norm(acc,2) > 8)
            [x,P] = mu_g(x, P, acc-mean_acc, var_acc, [0; 0; 9.82]);
            setAccDist(ownView,0);
            else
            setAccDist(ownView,1);
           end
    end

    mag = data(1, 8:10)';
    if ~any(isnan(mag))  % Mag measurements are available.
      % Do something

       abs(norm(mag,2) - norm(mean_mag,2));
          if (abs(norm(mag,2) - norm(mean_mag,2)) < 20)
            % earth magnetic field in world coordinates:
             m0 = [0 sqrt(mean_mag(1)^2 +mean_mag(2)^2) , mean_mag(3)]'; 
            [x,P] = mu_m(x, P, mag-mean_mag, var_mag, m0);
            setMagDist(ownView,0);
        else
            setMagDist(ownView,1);
        end
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
  figure(2)
  plot(q2euler(x)')
  hold on
  plot(q2euler(meas.orient)')
end
