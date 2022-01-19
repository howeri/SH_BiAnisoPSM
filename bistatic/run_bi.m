clc;clear;close all
%% aircraft example with 
% Bistatic RCS. Single freuqnecy 1GHz. 20 degree increment.
% Use bilinear equations to find SH coefficients.
% alpha^T * Qm * beta = rm
% Assume real gains
% Removed theta=0, and opposite directional data

%% Initial Setting
addpath('jmontalt-harmonicY-accdfe4')
M = 18*8*18*8; % # of data
K = 16;
f = 1e9;       % signal frequency
c = physconst('LightSpeed');
SH_degree = 12;

%% Import rcs data (Cut out extra data) (phi_obs, theta_obs, phi_exc, theta_exc)
rcs = importdata('../rcsData/bistatic_step20.txt');
rcs_Re = reshape(rcs.data(:,6),[18, 10, 18, 10]);  % phi directional RCS
rcs_Im = reshape(rcs.data(:,7),[18, 10, 18, 10]);  % phi directional RCS
obs_phi = reshape(rcs.data(:,4),[18, 10, 18, 10]); 
obs_theta = reshape(rcs.data(:,5),[18, 10, 18, 10]);
exc_phi = reshape(rcs.data(:,2),[18, 10, 18, 10]);
exc_theta = reshape(rcs.data(:,3),[18, 10, 18, 10]);

obs_phi = reshape(obs_phi(:, 2:9, :, 2:9), [18*8*18*8, 1]);
obs_theta = reshape(obs_theta(:, 2:9, :, 2:9), [18*8*18*8, 1]);
exc_phi = reshape(exc_phi(:, 2:9, :, 2:9), [18*8*18*8, 1]);
exc_theta = reshape(exc_theta(:, 2:9, :, 2:9), [18*8*18*8, 1]);
rcs_Re = reshape(rcs_Re(:, 2:9, :, 2:9)*2*sqrt(pi), [18*8*18*8, 1]); % 2*sqrt(pi) is for a CST offset 
rcs_Im = reshape(rcs_Im(:, 2:9, :, 2:9)*2*sqrt(pi), [18*8*18*8, 1]); % 2*sqrt(pi) is for a CST offset 

phi_data = 0:20:340;
theta_data = 20:20:160;

ind_del = [];
for p_out = phi_data
    for t_out = theta_data
        if wrapTo360(p_out+180)==360
            p_in_del = 0;
        else
            p_in_del = wrapTo360(p_out+180);
        end
        t_in_del = 180-t_out;       
        ind = find((obs_phi == p_out) & ...
                   (exc_phi == p_in_del) & ...
                   (obs_theta == t_out) & ...
                   (exc_theta == t_in_del));
        ind_del = [ind_del ind];
    end
end
ind_del = sort(ind_del);
ind_remain = 1:M;
ind_remain(ind_del)=[];
rcs_Re(ind_del)=inf;
rcs_Im(ind_del)=inf;
rcs_Re_mat = reshape(rcs_Re, [18, 8, 18, 8]);
rcs_Im_mat = reshape(rcs_Im, [18, 8, 18, 8]);
rcs_pow = rcs_Re.^2 + rcs_Im.^2;
rcs_pow_mat = reshape(rcs_pow, [18, 8, 18, 8]);

% figure
% subplot(1,2,1)
% plotImage(theta_data, phi_data, rcs_Re_mat(:,:,14,4), 'RCS from CST with \phi_{inc}=260^o & \theta_{inc}=80^oRCS (Re)')
% caxis('auto')
% subplot(1,2,2)
% plotImage(theta_data, phi_data, rcs_Im_mat(:,:,14,4), 'RCS from CST with \phi_{inc}=260^o & \theta_{inc}=80^oRCS (Im)')
% caxis('auto')
% figure
% plotImage(theta_data, phi_data, pow2db(rcs_pow_mat(:,:,14,4)), 'RCS power from CST  with \phi_{inc}=260^o & \theta_{inc}=80^o (dB)')
% caxis('auto')

%% Construct real spherical harmonics of degree SH_degree
SH_num = (SH_degree+1)^2; % SH_num=L
SH_matrix = [];
[PHI, THETA] = ndgrid(phi_data, theta_data);
for l=0:SH_degree
    for m=-l:l
        SH_matrix = cat(3, SH_matrix, harmonicY(l, m, deg2rad(THETA), deg2rad(PHI), 'type', 'real'));
        % Visualize spherical harmonics
%         figure
%         [x, y, z] = sph2cart(deg2rad(PHI), pi/2-deg2rad(THETA), abs(SH_matrix));
%         surf(x(:,:,end), y(:,:,end), z(:,:,end), SH_matrix(:,:,end));
%         axis equal
%         xlabel('X')
%         ylabel('Y')
%         zlabel('Z')
    end
end

%% Construct the unit direction matrix and combinations of adding 2 unit directions for delay computation
[PHI, THETA] = ndgrid(phi_data, theta_data);
rho = 1;
r_xy = rho .* sind(THETA);
x  = r_xy  .* cosd(PHI);
y  = r_xy  .* sind(PHI);
z  = rho .* cosd(THETA);
direction = cat(3, x, y, z);
direction2 = direction + reshape(permute(direction,[3,1,2]), [1 1 3 18 8]);
direction2 = permute(direction2, [1 2 4 5 3]);
direction2 = reshape(direction2, [18*8*18*8, 3]);

%% Construct E and Y matrixes of Q
ind_remain_complex = [ind_remain M+ind_remain]; 
addpath('../monostatic/results')
load('sc16SH12_locations.mat')
locations = [X_sc Y_sc Z_sc];
lead = locations*direction2'/c;
e_Re = real(exp(1i*2*pi*f*lead));
e_Im = imag(exp(1i*2*pi*f*lead));
e = cat(2, e_Re, e_Im);
e = e(:,ind_remain_complex);

SH_temp = reshape(SH_matrix, [18*8, SH_num]);
SHob_temp = repmat(SH_temp, 18*8, 1);
SHob_temp2 = reshape(SHob_temp', [1, SH_num, M]);
Y_obs_temp = repmat(SHob_temp2, SH_num, 1, 1);
Y_obs_temp2 = cat(3, Y_obs_temp, Y_obs_temp);
Y_obs = Y_obs_temp2(:,:,ind_remain_complex);

SHex_temp = repelem(SH_temp, 18*8, 1);
SHex_temp2 = reshape(SHex_temp', [SH_num, 1, M]);
Y_exc_temp = repmat(SHex_temp2, 1, SH_num, 1);
Y_exc_temp2 = cat(3, Y_exc_temp, Y_exc_temp);
Y_exc = Y_exc_temp2(:,:,ind_remain_complex);

%% Solve the Bilinear Problem
r = [rcs_Re(ind_remain); rcs_Im(ind_remain)];
tic;
[a, b] = NorIterAlg(e, Y_obs, Y_exc, r, K);
toc
%% Plot Modeled RCS
% addpath('results')
% load('sc16SH12_NorIterAlg_sol.mat')
rcs_Re_modeled_plot = inf(18*8, 1);
rcs_Im_modeled_plot = inf(18*8, 1);
ind = reshape(1:(18*8*18*8), [18,8,18,8]);
ind = reshape(ind(:,:,14,4), [18*8, 1]); % select interested area
e = exp(1i*2*pi*f*lead);
for i=1:(18*8)   
    if ~ismember(ind(i),ind_del)            
        SH_temp = reshape(SH_matrix, [18*8, SH_num]);       
        SHob_temp = repmat(SH_temp, 18*8, 1);
        SHob_temp2 = SHob_temp(ind(i), :);
        Y_obs = repmat(SHob_temp2, SH_num, 1);
        SHex_temp = repelem(SH_temp, 18*8, 1);
        SHex_temp2 = SHex_temp(ind(i), :);
        Y_exc = repmat(SHex_temp2', 1, SH_num);
        Q_Re_plot = [];
        Q_Im_plot = [];
        for k=1:K
            e_k_plot = e(k,ind(i));
            Q_Re_plot = blkdiag(Q_Re_plot, Y_exc.*Y_obs.*real(e_k_plot));
            Q_Im_plot = blkdiag(Q_Im_plot, Y_exc.*Y_obs.*imag(e_k_plot));
        end          
        rcs_Re_modeled_plot(i) = a'*Q_Re_plot*b;
        rcs_Im_modeled_plot(i) = a'*Q_Im_plot*b;
    end
end
rcs_Re_modeled_mat_plot = reshape(rcs_Re_modeled_plot, [18, 8]);
rcs_Im_modeled_mat_plot = reshape(rcs_Im_modeled_plot, [18, 8]);

figure
subplot(1,2,1)
plotImage(theta_data, phi_data, rcs_Re_modeled_mat_plot, 'Modeled RCS with \phi_{inc}=260^o & \theta_{inc}=80^oRCS (Re)')
subplot(1,2,2)
plotImage(theta_data, phi_data, rcs_Im_modeled_mat_plot, 'Modeled RCS with \phi_{inc}=260^o & \theta_{inc}=80^oRCS (Im)')
figure
plotImage(theta_data, phi_data, pow2db(rcs_Re_modeled_mat_plot.^2 + rcs_Im_modeled_mat_plot.^2), 'Modeled RCS power with \phi_{inc}=30^o & \theta_{inc}=60^o (dB)')

%% Plot Point Scatterer Model
SH_num = (SH_degree+1)^2; % SH_num=L
SH_matrix = [];
[PHI, THETA] = ndgrid(1:360, 1:180);
for l=0:SH_degree
    for m=-l:l
        SH_matrix = cat(3, SH_matrix, harmonicY(l, m, deg2rad(THETA), deg2rad(PHI), 'type', 'real'));
    end
end

figure
X = sind(THETA).*cosd(PHI);
Y = sind(THETA).*sind(PHI);
Z = cosd(THETA);
a_mat = reshape(a, [SH_num, K]);
for k=1:K
    X_k = X + X_sc(k);
    Y_k = Y + Y_sc(k);
    Z_k = Z + Z_sc(k);   
    a_k = reshape(a_mat(:,k), [1,1,SH_num]);
    surf(X_k, Y_k, Z_k, (sum(SH_matrix .* a_k, 3)))
    hold on
end
colorbar
title('Point Scatterer Model (incident angle)')
caxis([0 0.2])
set(gca,'FontSize',20)
shading interp
axis square

figure
X = sind(THETA).*cosd(PHI);
Y = sind(THETA).*sind(PHI);
Z = cosd(THETA);
b_mat = reshape(b, [SH_num, K]);
for k=1:K
    X_k = X + X_sc(k);
    Y_k = Y + Y_sc(k);
    Z_k = Z + Z_sc(k);   
    b_k = reshape(b_mat(:,k), [1,1,SH_num]);
    surf(X_k, Y_k, Z_k, (sum(SH_matrix .* b_k, 3)))
    hold on
end
colorbar
title('Point Scatterer Model (scatter angle)')
caxis([0 2000])
set(gca,'FontSize',20)
shading interp
axis square

%% Plot cost while varying L
r1 = [187.3806, 182.5388, 178.5335, 172.2588, 165.3059, 159.0383, 151.4727, 146.4947, 142.1175, 137.2593];
r2 = [183.9354, 177.8689, 171.0827, 161.8061,154.3891, 144.0823, 136.8742, 130.1228, 123.9262, 119.3853];
figure
hold on
plot(3:12, r1, '-o', 'LineWidth', 2)
plot(3:12, r2, '-o', 'LineWidth', 2)
set(gca,'FontSize',20)
legend('with 16 scatterers', 'with 24 scatterers')
xlabel('Spherical Harmonic Degree (L)')
ylabel('Cost (J)')
