    % ========================================================================
% This script contains the calibration parameters of the setup.
% ========================================================================

%% Coils 
% Dor B (ALL NOT CALIBRATED YET!):
% Bx, By, NMR: 0.613 G/A, with CS @ 10 mA/V, calibration from 19/11/23
calib.coils.dBy_dI = 0.613; % G/A
calib.coils.dBx_dI = 0.613; % G/A
calib.coils.dBNMR_dI = 0.613; % G/A
calib.coils.dIy_dV = 10e-3; %A/V
calib.coils.dIx_dV = 10e-3; %A/V
calib.coils.dINMR_dV = 0.1e-3; %A/V
calib.coils.dBx_dVx = calib.coils.dBx_dI * calib.coils.dIx_dV ; % G/V
calib.coils.dBy_dVy = calib.coils.dBy_dI * calib.coils.dIy_dV;% 6.13e-3   ; % G/V
calib.coils.dBnmr_dVnmr = calib.coils.dBNMR_dI * calib.coils.dINMR_dV; % G/V

% Bz: 2.98 G/A, calibration from 19/11/23
calib.coils.dBz_dIz = 2.98;

% ESR: 1.9 G/A, calibration from 19/11/23
calib.coils.dBesr_dVesr = 11; % G/V

% WP.coils.dBz_dVz_feedbackCoil = (578+588)/2*1e-6; % Roy's calibration from 11/04/21 & 17/04/21
% Mini-Me:
% WP.coils.dBx_dVx = 25.74e-3; % G/V
% WP.coils.dBy_dVy = 28.69e-3; % G/V [Roy Elkabetz - 25/05/2022]
% WP.coils.dBnmr_dVnmr = 1.241e-3; % G/V
% WP.coils.dBz_dIz = 10.76;
% WP.coils.dBesr_dVesr = 11; % G/V
% WP.coils.dBz_dVz_feedbackCoil = (578+588)/2*1e-6; % Roy's calibration from 11/04/21 & 17/04/21
% WP.coils.dgradBzI_dV = 2; % mA/V


%%
prm.calib = calib;
clear calib