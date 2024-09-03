% Compute TBG local density of states using KPM with Jackson smoothing. 
% Local Chebyshev weights <v|T_n(H)|v> computed using get_cheb_wgts_ldos script.
% We load a vector cheb_wgts of weights produced by that script.

% Input parameters
filename = 'r100_p1000_ldos.mat';
p = 1000;     % Chebyshev degree
m = 6;        % Order of method with respect to broadening parameter eta
eta = 0.03;   % Broadening parameter
dE = 0.01;    % Energy grid spacing

addpath('hodc','hodc/kernels');

load(['cheb_wgts_data/',filename]); % Load parameters and Chebyshev weights from file

E = (-E_range):dE:E_range; % Energy grid

% Compute local densities of states
ldos = hodc_ldos(m, eta, p, E/E_range, cheb_wgts(1:p));

%addpath ~/Documents/MATLAB/export_fig/

figure(1);
plot(E, ldos, '.-');
%export_fig('dos_hodc_1000.pdf');