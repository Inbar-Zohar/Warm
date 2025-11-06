% ========================================================================
%   This script is measuring Alkali FID and computing the T2 and resonance 
%   frequncies from the data. The code assumes Bz is set, Bx and By are 
%   zeroed, BPD is balanced by a lambda plate. 
%   The signal from BPD is connected ch1 scope1. 
% ========================================================================
%
%   TODO:
%   # compute the larmour frequency from known values instead of reading it
%   from the ESR function generator
%   # set scope to averaging so that it trnafers 12bit data (check it)
%   #set Bx to DC
%
%   Set scope scale dynamically - Roy 23/09/2021
%   Corrected the fitting routine - CF 29/09/21
%
% ========================================================================


% Disable PID 129 131 ESR

% Disable NMR 129 and 131 and ESR
AG1.OutputOFF(1);   AG2.OutputOFF(1);   
AG3.OutputOFF(1); 
AG6.OutputOFF(1);
% Disable Bx
% AG1.OutputOFF(2); 

% Set  AG6 ch2 to trigger all with a square pulse
AG6.Square(2, 133, 1.5, 1.5, 0); AG6.Termination(2,50);  AG6.OutputON(2);

% read larmour frequency of Alkali from function generator 
% [stt, wvfm] = AG3.Output(1, [], true);
% prm.fresESR = wvfm{1}.p(1);
prm.fresESR = 150e3; %Manually set (temporarily)

% % compute the larmour frequency from known values instead of reading the
% % ESR function generator
% prm.fresESR = prm.WP.coils.dBz_dIz * prm.WP.fields.Iz * c.gRb85;

% Set AG2 ch2 to By pulses, with ext trig
dt = 1 / 10 / prm.fresESR; % dt is shorter 10 times than larmours cycle. (NMR limit)
pulse_amp = 3.5;

% warning('setting By pulse parameters to be OCFES like');
% dt = 5e-6;
% pulse_amp = 0.67;
 
DorB.Bx.apply('DC',prm.WP.fields.Bx.B0);
AG2.BurstPulse(2, 1e3, pulse_amp, dt, prm.WP.fields.By.B0/10-pulse_amp/2, 0);

% Set scope 1 ch1 with trigger on ch1
% Configure scope to take data
chnl = 1;
Tscale = 20e-6; % sec per dev
Vscale = 1; % volt per dev

% set scope
curr_inst = scope2;
% curr_inst.LoadDefaultSetup;
curr_inst.Stop;

% Coupling
curr_inst.setChCoupling(chnl, 'DC');
curr_inst.setChImpedance(chnl, 'M'); % 50 Ohm (F), or 1 MOhm (M)

% Acquire
% curr_inst.HighRes;
% curr_inst.AverageAQN(2); % for 12bit data (not working as is)

% horizontal
curr_inst.setTscale(Tscale);
curr_inst.setTdelay(4.00 * Tscale);

% vertical offset
if 1
    curr_inst.TrigMode('NORM');
    curr_inst.TrigSource('EXT');
    curr_inst.TrigSlope('POS');
    curr_inst.TrigThresh(1.0);
    
    curr_inst.autoVrange(chnl,Tscale,[],0);
    pause(0.5)
    curr_inst.autoVrange(chnl,Tscale,[],0);
    curr_inst.Single;
    curr_inst.readyToRead(10);   
    [t, v, ~] = curr_inst.Read(chnl);
    
else
    curr_inst.setVoffset(chnl, 0);
    
    % setting the scope vscale dynamically
    try_counter = 0;
    in_once_flag = 0;
    out_once_flag = 0;
    
    while try_counter < 30 && ~in_once_flag && Vscale >= 0.002
        % setting signal to scope bounds by enlarging the scale until signal is
        % out of scope and than setting it to smaller scale until signal is in
        % again.
        
        % vertical scale
        curr_inst.setVscale(chnl, Vscale);
        curr_inst.Run;
        
        % Trigger
        curr_inst.TrigMode('NORM');
        curr_inst.TrigSource('EXT');
        curr_inst.TrigSlope('POS');
        curr_inst.TrigThresh(1.0);
        curr_inst.Single;
        
        % read data from scope
        curr_inst.readyToRead(10);
        [t, v, ~] = curr_inst.Read(chnl);
        
        in_flag = curr_inst.sigInbounds(v, chnl);
        if ~in_flag
            % signal not contained in scope, enlarge scale
            Vscale = Vscale * 1.25;
            try_counter = try_counter + 1;
            out_once_flag = 1;
            
        elseif in_flag && ~out_once_flag
            % signal is contained in scope from the start, set smaller scale
            Vscale = Vscale * 0.5;
            try_counter = try_counter + 1;
            
        elseif in_flag && out_once_flag
            in_once_flag = 1;
        end
        
        
    end
end

% Release scope from Single mode
curr_inst.Run;

meas.FID_Rb_t = t;
meas.FID_Rb_v = v;

% % uncomment to plot waveform 
% meas=create_fig(meas); plot(meas.FID_Rb_t / Tscale, meas.FID_Rb_v); xlabel('time'); ylabel('amplitude') 

clear t v 

%% FID rubidium Analyze


% logger.log('Analyzing measurement results')
Alkali_FID_measurement_analysis;
