clc;clear;close all
%% Aircraft Example 
% Monostatic RCS. Single freuqnecy 1GHz. 10 degree increment.
% Use linear programming with sparsity constraints to find SH coefficients.
% S * alpha = r
% Assume real gains

%% Initial Setting
addpath('jmontalt-harmonicY-accdfe4')
M = 18*36; % # of data
K = 10*10*10; % number of possible scatterer locations
f = 1e9;   % signal frequency
c = physconst('LightSpeed');
SH_degree = 12;
N_scatterers = 16;

%% Import rcs Data (Cut out extra data)
rcs = importdata('../rcsData/monostatic_step10.txt');
phi_data = rcs.data(1:M,2);
theta_data = rcs.data(1:M,3);
rcs_Re = rcs.data(1:M,6)*2*sqrt(pi); % 2*sqrt(pi) is from a CST offset 
rcs_Im = rcs.data(1:M,7)*2*sqrt(pi); % 2*sqrt(pi) is from a CST offset 
rcs_pow = reshape(abs(rcs_Re+1i*rcs_Im).^2, [36,18]);
figure
plotImage(0:10:170, 0:10:350, rcs_pow, 'RCS from CST')
figure
plotImage(0:10:170, 0:10:350, pow2db(rcs_pow), 'RCS from CST (dB)')

%% Construct Real Spherical Harmonics
SH_num = (SH_degree+1)^2; % SH_num=L
SH_matrix = [];
[PHI, THETA] = ndgrid(0:10:350, 0:10:170);
for l=0:SH_degree
    for m=-l:l
        SH_matrix = cat(3, SH_matrix, harmonicY(l, m, deg2rad(THETA), deg2rad(PHI), 'type', 'real'));
%         % Visualize spherical harmonics
%         figure
%         [x, y, z] = sph2cart(deg2rad(PHI), pi/2-deg2rad(THETA), abs(SH_matrix));
%         surf(x(:,:,end), y(:,:,end), z(:,:,end), SH_matrix(:,:,end));
%         axis equal
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
    end
end
SH = reshape(SH_matrix, [36*18, SH_num]);

%% Construct a list of direction unit vectors with theta_data and phi_data
rho = 1;
r_xy = rho .* sind(theta_data);
x_dir  = r_xy  .* cosd(phi_data);
y_dir  = r_xy  .* sind(phi_data);
z_dir  = rho .* cosd(theta_data);
direction = [x_dir y_dir z_dir]';

%% Construct S Matrix
[X_sc, Y_sc, Z_sc] = ndgrid(1:10 ,1:10, 1:10);  %possible scatterer locations
X_sc = reshape(X_sc, [1000 1]);
Y_sc = reshape(Y_sc, [1000 1]);
Z_sc = reshape(Z_sc, [1000 1]);
locations = [X_sc Y_sc Z_sc];
lead = 2*locations*direction/c;
S_temp = transpose(exp(1i*2*pi*f*lead));
S = repelem(S_temp, 1, SH_num) .* repmat(SH, 1, 1000);
S_stacked = [real(S); imag(S)];

%% Solve Sparse Linear Problem
tic
[X_hat, b_hat, res] = OMP_modified(S_stacked, [rcs_Re; rcs_Im], N_scatterers, SH_num); 
toc
b_hat = b_hat(1:M)+ 1i*b_hat(M+1:end);
disp(['Residual: ',num2str(res)]); % norm of the residual vector
RMSE_dB = sqrt (sum((pow2db(abs(rcs_Re+1i*rcs_Im).^2) - pow2db(abs(b_hat).^2)).^2) / length(b_hat));
disp(['RMSE of powers in dB scale: ',num2str(RMSE_dB)]);

%% Plot Comparisons
figure
subplot(2,2,1)
plotImage(0:10:170, 0:10:350, reshape(rcs_Re, [36,18]), 'RCS from CST (Re)')
caxis([-200 100])
subplot(2,2,2)
plotImage(0:10:170, 0:10:350, reshape(rcs_Im, [36,18]), 'RCS from CST (Im)')
caxis([-300 100])
subplot(2,2,3)
plotImage(0:10:170, 0:10:350, reshape(real(b_hat), [36,18]), 'Modeled RCS (Re)')
caxis([-200 100])
subplot(2,2,4)
plotImage(0:10:170, 0:10:350, reshape(imag(b_hat), [36,18]), 'Modeled RCS (Im)')
caxis([-300 100])

figure
subplot(1,2,1)
plotImage(0:10:170, 0:10:350, pow2db(rcs_pow), 'RCS from CST (dB)')
caxis([-20 50])
subplot(1,2,2)
plotImage(0:10:170, 0:10:350, pow2db(reshape(abs(b_hat).^2, [36,18])), 'Modeled RCS (dB)')
caxis([-20 50])

%% Plot Point Scatterer Model
figure
X = sind(THETA).*cosd(PHI);
Y = sind(THETA).*sind(PHI);
Z = cosd(THETA);
X_hat_mat = reshape(X_hat, [SH_num, 10, 10, 10]);
sol_extract = squeeze(sum(X_hat_mat, 1));
[X_sc,Y_sc,Z_sc] = ind2sub([10,10,10],find(sol_extract));
reflectionGain_max = [];
reflectionGain_min = [];
for k=1:N_scatterers
    X_k = X + X_sc(k);
    Y_k = Y + Y_sc(k);
    Z_k = Z + Z_sc(k);
    X_hat_mat_ = reshape(X_hat_mat(:, X_sc(k), Y_sc(k), Z_sc(k)), [1,1,SH_num]);
    reflectionGain_max = [reflectionGain_max max(abs(sum(SH_matrix .* X_hat_mat_, 3)),[],'all')];
    reflectionGain_min = [reflectionGain_min min(abs(sum(SH_matrix .* X_hat_mat_, 3)),[],'all')];
    surf(X_k, Y_k, Z_k, abs(sum(SH_matrix .* X_hat_mat_, 3)))
    hold on
end
axis equal
colorbar
disp([ 'Max abs of reflection gains: ',num2str(max(reflectionGain_max))]);
disp([ 'Max abs of reflection gains: ',num2str(min(reflectionGain_min))]);

%% Save Locations
X_hat_mat = reshape(X_hat, [SH_num, 10, 10, 10]);
sol_extract = squeeze(sum(X_hat_mat, 1));
[X_sc, Y_sc, Z_sc] = ind2sub([10,10,10], find(sol_extract));
% save('sc16SH12_locations','X_sc','Y_sc','Z_sc');
