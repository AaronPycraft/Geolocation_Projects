% Project
clear; clc;
format shortEng; format compact;
%%
c = 3e8;
%%

rx1 = [39.172739,-76.897929];
rx2 = [39.171851,-76.898607];
tx  = [39.167385,-76.899579];
rx3 = [39.165859,-76.899764];
rx4 = [39.164740,-76.900857];
rx5 = [39.161708,-76.901038];
loc = [rx1 rx2 tx rx3 rx4 rx5];

tx_indx = 3;

range(1) = ll2km(rx1, tx);
range(2) = ll2km(rx2, tx);
range(3) = ll2km(tx, tx);
range(4) = ll2km(rx3, tx);
range(5) = ll2km(rx4, tx);
range(6) = ll2km(rx5, tx);

% range = 650;
f = 430e6;
lambda = c/f;
FSPL_dB = 20*log10(4*pi*f/c*range)

dt = range/c

toa = dt - dt(1)
%%
% constants
kB=physconst('Boltzmann');
T0 = 290;
f = 430e6;
lambda = c/f;
% Tx
Pt = 7;
Gt = 3;
% Rx
Gr = 0;
NFr = 3.5;
Br = 4e6;
% Calculate
MDS = 10*log10(kB*T0/(1e-3)) + NFr + 10*log10(Br)
MDS = -100

Rmax = sqrt(10^((Gt+Gr+Pt-MDS)/10)) * lambda/(4*pi)



