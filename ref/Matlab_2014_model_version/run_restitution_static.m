function [APD_restitution_static] = run_restitution_static(modelswitch)
%**************************************************************************
% Static restitution simulation protocol
%  > starting from a pacing steady-state (pacing_1AP_ss_BCL_1000.dat)
%  > increasing BCL = 250:10:1000;
%  > 1 AP simulations at each BCL for calculating the results
%--------------------------------------------------------------------------
% 04.10.2012, Jussi Koivumäki
%**************************************************************************

%% Stimulus for s2
sim_length = 1;
BCL_s2 = 1000;
stim = stim_creator_BCL(BCL_s2, 1e-3, sim_length);

%% Load s1
s1 = load('pacing_1AP_ss_BCL_1000.dat');
time = s1(:, 1);
V = s1(:, 2);

%% Calculate time point of APD90 in s1
V_max = max(V);
index_AP_huippu = find(V == V_max, 1);
time_shifted = time(index_AP_huippu:end);
V_shifted = V(index_AP_huippu:end);
V_min = min(V_shifted);
V_ampl = V_max - V_min;
time_90 = time_shifted(find(V_shifted <= (V_max - 0.90*V_ampl), 1));

%% Calculate APD restitution
if strcmp(modelswitch, 'nSR')
%    BCL_sequence = 1000 - (0:10:750); % BCL = [1000, 250] ms
    BCL_sequence = 1000 - (0:10:850); % BCL = [1000, 150] ms
else
%    BCL_sequence = 1000 - (0:10:780); % BCL = [1000, 220] ms
    BCL_sequence = 1000 - (0:10:850); % BCL = [1000, 150] ms
end
DI = zeros(1, length(BCL_sequence));
APD_90 = zeros(1, length(BCL_sequence));
Nass = zeros(1, length(BCL_sequence));
CaTdias = zeros(1, length(BCL_sequence));
NCX_fwd = zeros(1, length(BCL_sequence));
NCX_rev = zeros(1, length(BCL_sequence));

Vcytosol = [0.506253 1.51876 2.53127 3.54377]; % Without ss compartment
w = Vcytosol'/sum(Vcytosol) * 1000; % volumetric

for k = 1 : length(BCL_sequence);
    DI(k) = BCL_sequence(k) - time_90*1000;
    % Define initial values for s2
    initcond_index = find(time >= BCL_sequence(k)/1000, 1);
    fid = fopen('y0_nosave.dat', 'w');
    style = '';
    for n = 1:(29 + 6 + 4*2) % functions, differential variables and time
        style = [style ' %6.6e'];
    end
    fprintf(fid, style, s1(initcond_index, 2:44));
    fclose(fid);
    % Run 1 AP simulations
    main_human_atrial(stim, 'y0_nosave.dat', 'nosave.dat', sim_length, BCL_s2, 0, modelswitch);
    s2 = load('pacing_nosave.dat');
    time = s2(:, 1);
    V = s2(:, 2);
    % Calculate APD90 in s2
    V_max = max(V);
    index_AP_huippu = find(V == V_max, 1);
    time_AP_huippu = time(find(V == max(V), 1));
    time_shifted = time(index_AP_huippu:end);
    V_shifted = V(index_AP_huippu:end);
    V_min = min(V_shifted);
    V_ampl = V_max - V_min;
    APD_90(k) = time_shifted(find(V_shifted <= (V_max - 0.90*V_ampl), 1)) - time_AP_huippu;
    Nass(k) = mean(s2(:, 33));
    CaT = s2(:, 37)*w(1) + s2(:, 38)*w(2) + s2(:, 39)*w(3) + s2(:, 40)*w(4);
    CaTdias(k) = min(CaT);
    NCX_fwd(k) = min(s2(:, 56));
    NCX_rev(k) = max(s2(:, 56));
end

APD_restitution_static = [DI' (APD_90*1000)'];
save 'APD_restitution_static.mat' APD_restitution_static;

Nass_static = [BCL_sequence' Nass']; save 'Nass_static.mat' Nass_static;
APD90_static = [BCL_sequence' (APD_90*1000)']; save 'APD90_static.mat' APD90_static;
CaTdias_dynamic_static = [BCL_sequence' CaTdias']; save 'CaTdias_dynamic_static.mat' CaTdias_dynamic_static;
NCX_dynamic_static = [BCL_sequence' NCX_fwd' NCX_rev']; save 'NCX_dynamic_static.mat' NCX_dynamic_static;

%% Plot APD restitution
%figure;
%plot(APD_restitution_static(:, 1), APD_restitution_static(:, 2), 'x');
%xlabel('DI (ms)');
%ylabel('APD_{90} (ms)');
%title('Static APD restitution');

end