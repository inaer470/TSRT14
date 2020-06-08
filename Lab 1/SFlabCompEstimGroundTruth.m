function SFlabCompEstimGroundTruth(trajEst,micPos)
%
% Function that plots the estimated trajectory and microphone positions on
% top of an image of the ground truth.
%
% function SFlabCompEstimGroundTruth(trajEst,micPos)
%
% INPUTS:
% trajEst   -  2XN-vector with the estimated trajectory, first row is
%              X-coordinate, second row is Y-coordinate.
% micPos    -  2xN-matrxi with the microphone positions, first row is
%              X-coordinate, second row is Y-coordinate.
%
%

% Martin Skoglund and Karl Granström
% 2009-03-24

% Rotation matrix
R =[...
   0.999774521898456  -0.000924915849553   0.021214379401393;...
   0.001705399360543  -0.992326335159808  -0.123634688341604;...
   0.021165939046876   0.123642990415857  -0.992101009950745];
% Translation matrix
T =[...
  -0.583118758321543;...
   0.458078043186633;...
   1.412654371966365];
% Focal lengths
f =1.0e+003*[...
   2.389065133857016;...
   2.393358741136121];
% Principal point coordinates
c =1.0e+003*[...
   1.120043265067930;...
   0.861829343480461];

% coordinates in camera frame
g = [R T; 0 0 0 1];
PI0 = [eye(3) zeros(3,1)];
xc =PI0*g*[trajEst; zeros(1,size(trajEst,2)) ; ones(1,size(trajEst,2))];
mc =PI0*g*[micPos ; zeros(1,size(micPos,2)) ; ones(1,size(micPos,2))];

% normalized coordinates
xn(1,:) = xc(1,:)./xc(3,:);
xn(2,:) = xc(2,:)./xc(3,:);
mn(1,:) = mc(1,:)./mc(3,:);
mn(2,:) = mc(2,:)./mc(3,:);

% pixel coordinates
K = [
    f(1)  0     c(1);
    0     f(2)  c(2);
    0     0     1];
xp = K*[xn ; ones(1,size(trajEst,2))];
mp = K*[mn ; ones(1,size(micPos,2))];

% Load image of ground truth
I=imread('SFlabGroundTruth.jpg');
I = rgb2gray(I);

% Compute axis min and max
[Iy,Ix] = size(I);
extra = 50;
xmin = min([0 xp(1,:) mp(1,:)])-extra;
xmax = max([Ix xp(1,:) mp(1,:)])+extra;
ymin = min([0 xp(2,:) mp(2,:)])-extra;
ymax = max([Iy xp(2,:) mp(2,:)])+extra;

% Plot
figure
imagesc(I)
colormap(gray)
hold on
plot(xp(1,:),xp(2,:),'r.','markersize',6,'linewidth',2)
plot(mp(1,:),mp(2,:),'b.','markersize',6,'linewidth',2)
hold off
title('Estimated trajectory in red, microphone positions in blue')
axis image
axis([xmin xmax ymin ymax])