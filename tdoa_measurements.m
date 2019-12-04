% HW1 - Aaron Pycraft
clear all; close all;
format compact; format shortEng;
clc;
%% Constants
us = 1e-6;
km = 1e3;
c = 299792.458; % km/s
%% Inputs
% tx = [39.167385,-76.899579];
% 
% loc= [39.172739,-76.897929;...
%       39.171851,-76.898607;...
%       39.165859,-76.899764;...
%      39.161708,-76.901038];
% %{
%Measurements
tx =  [ -76.899704, 39.167397];
loc = [ -76.900037, 39.167113;...
        -76.899762, 39.167898;...
        -76.898859, 39.167249;...
        -76.899312, 39.167308];
% %}
% convert lat long to xy grid projection
scaleFactor = 27.4786e6;
% scaleFactor(1) = 1;%36.4206/1.1845e-006;
% scaleFactor(2) = 1;%36.4206/1.1845e-006; %* 131.0325e+000/156.0120e+000;
% scaleFactor(3) = 1;%36.4206/1.1845e-006;
% scaleFactor(4) = 1;%36.4206/1.1845e-006;
[locXs(1), locYs(1)] = Spherical2AzimuthalEquidistant(loc(1,2), loc(1,1), tx(2), tx(1), 0, 0, scaleFactor);
[locXs(2), locYs(2)] = Spherical2AzimuthalEquidistant(loc(2,2), loc(2,1), tx(2), tx(1), 0, 0, scaleFactor);
[locXs(3), locYs(3)] = Spherical2AzimuthalEquidistant(loc(3,2), loc(3,1), tx(2), tx(1), 0, 0, scaleFactor);
[locXs(4), locYs(4)] = Spherical2AzimuthalEquidistant(loc(4,2), loc(4,1), tx(2), tx(1), 0, 0, scaleFactor);


%{
% Theoretical
loc = [ 39.171828,-76.897896;...
        39.168604,-76.899255;...
        39.167339,-76.899604;...
        39.164753,-76.900531;...
        39.161945, -76.900950];
tx = [39.167339,-76.899604];
tar_indx = 3;
        %}


for ii=1:size(loc,1)
%     range(ii) = ll2km([loc(ii,:)], [loc(tar_indx,:)]);
    rangeFromTX(ii) = ll2km([loc(ii,:)], [tx]); % lat-long to meters
    rangeFromLoc1(ii) = ll2km([loc(ii,:)], loc(1,:)); % lat-long to meters
end

% Simulated TOA
toa = rangeFromTX/c;

% measured TOA



%% example
%{
P1 = [-1.5  -2.0];
P2 = [ 2.0   2.0];
P3 = [-2.5   2.5];
P4 = [ 2.0  -1.0];

td21 =  0.095 * us
td31 =  4.950 * us
td41 = -1.880 * us
%}
%%
P1 = loc(1,:);
P2 = loc(2,:);
P3 = loc(3,:);
P4 = loc(4,:);
% P5 = loc(5,:);
tdoa = toa-toa(1);
td21 = toa(2)-toa(1);
td31 = toa(3)-toa(1);
td41 = toa(4)-toa(1);


%% Calculations
syms xi yi xj yj xt yt c T;
% RHS=sqrt((xt-xi)^2+(yt-yi)^2) - sqrt((xt-xj)^2+(yt-yj)^2);
% soln = solve(RHS-c*T, yt)
soln =  [(xi^2*yj - xi^2*yi + xj^2*yi - xj^2*yj + yi*yj^2 + yi^2*yj - yi^3 - yj^3 + T*c*((- T^2*c^2 + xi^2 - 2*xi*xj + xj^2 + yi^2 - 2*yi*yj + yj^2)*(- T^2*c^2 + xi^2 + 2*xi*xj - 4*xi*xt + xj^2 - 4*xj*xt + 4*xt^2 + yi^2 - 2*yi*yj + yj^2))^(1/2) + T^2*c^2*yi + T^2*c^2*yj + 2*xi*xt*yi - 2*xi*xt*yj - 2*xj*xt*yi + 2*xj*xt*yj)/(2*T^2*c^2 - 2*yi^2 + 4*yi*yj - 2*yj^2);...
         (xi^2*yj - xi^2*yi + xj^2*yi - xj^2*yj + yi*yj^2 + yi^2*yj - yi^3 - yj^3 - T*c*((- T^2*c^2 + xi^2 - 2*xi*xj + xj^2 + yi^2 - 2*yi*yj + yj^2)*(- T^2*c^2 + xi^2 + 2*xi*xj - 4*xi*xt + xj^2 - 4*xj*xt + 4*xt^2 + yi^2 - 2*yi*yj + yj^2))^(1/2) + T^2*c^2*yi + T^2*c^2*yj + 2*xi*xt*yi - 2*xi*xt*yj - 2*xj*xt*yi + 2*xj*xt*yj)/(2*T^2*c^2 - 2*yi^2 + 4*yi*yj - 2*yj^2)];
% Plug in values
xt = [-200:1:200];
c = 299792.458; % km/s
% Convert system to xy grid

xref = 0; yref = 0;  % coordinate system origin
x1 = locXs(1); y1 = locYs(1);
x2 = locXs(2); y2 = locYs(2);
x3 = locXs(3); y3 = locYs(3);
x4 = locXs(4); y4 = locYs(4);

% tdoa reference loc
xj = x1; yj = y1;

% Hyperbole of P1, P2
T = td21;
xi = x2;
yi = y2;
result12 = eval(soln);
hyp12_1 = result12(1,:);
hyp12_2 = result12(2,:);

% Hyperbole of P3, P1
T = td31;
xi = x3;
yi = y3;
result31 = eval(soln);
hyp13_1 = result31(1,:);
hyp13_2 = result31(2,:);

% Hyperbole of P1, P4
T = td41;
xi = x4;
yi = y4;
result41 = eval(soln);
hyp14_1 = result41(1,:);
hyp14_2 = result41(2,:);


%% Plots
figure(1); 
grid minor; 
hold on;
axis tight;
%% plot hyperboles

% plot(xt, hyp12_1, '--c');
plot(xt, hyp12_2, '-.c');

plot(xt, hyp13_1, '--b');
% plot(xt, hyp13_2, '-.b');

% plot(xt, hyp14_1, '--g');
plot(xt, hyp14_2, '-.g');
legend('hyperbole 12', 'hyperbole 13', 'hyperbole 14');

%% Plot Base Stations
plot(xref, yref, 'bx', 'HandleVisibility','off'); 
plot(x1, y1, 'ro', 'HandleVisibility','off'); 
plot(x2, y2, 'ro', 'HandleVisibility','off'); 
plot(x3, y3, 'ro', 'HandleVisibility','off'); 
plot(x4, y4, 'ro', 'HandleVisibility','off'); 

zoom = 30;
xlim([min(locXs)-zoom , max(locXs)+zoom ]);
ylim([min(locYs)-zoom , max(locYs)+zoom ]);

text(xref+0.001e3, yref, 'TX');
text(x1+0.001e3, y1, 'P1');
text(x2+0.001e3, y2, 'P2');
text(x3+0.001e3, y3, 'P3');
text(x4+0.001e3, y4, 'P4');

title('Scenario Map');
xlabel('Longitude coordinate (m)');
ylabel('Latitude coordinate (m)');

%{
% plot(tx(1), tx(2), 'bx'); 
% plot(x1, y1, 'ro'); 
% plot(x2, y2, 'ro'); 
% plot(x3, y3, 'ro'); 
% plot(x4, y4, 'ro'); 

% xlim([min(loc(:,1))-0.5e-3, max(loc(:,1))+0.5e-3]);
% ylim([min(loc(:,2))-0.5e-3, max(loc(:,2))+0.5e-3]);

% text(tx(1)+150e-6, tx(2), 'TX');
% text(x1+100e-6, y1, 'P1');
% text(x2+100e-6, y2, 'P2');
% text(x3+100e-6, y3, 'P3');
% text(x4+100e-6, y4, 'P4');

% title('Scenario Map');
% xlabel('Longitude coordinate (deg)');
% ylabel('Latitude coordinate (deg)');
%}
% set(gca, 'fontsize', 16);

%% Error Ellipse
initEst = [-4, -1.6854];
P = cov([xref, yref].', initEst.') + cov(initEst, [xref, yref])
ellipse = construct_ellipse(1000, initEst', P, chi2inv(0.95, 2));
plot(ellipse(1,:), ellipse(2,:))