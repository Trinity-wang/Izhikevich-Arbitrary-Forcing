%% ARBITRARY FORCING OF IZHIKEVICH SPIKING NEURON MODEL 
%  COPYWRITE: ANDREW HANSEN, 2015
%  BRANDEIS UNIVERSITY

%% 1. SIMULATION PARAMETERS (PLAY AROUND WITH THESE)

clc;
clear all;
close all;

sim_bnd = [0, 1]; % Initialize simulation time interval (seconds)
sim_res = 0.001; % Initialize simulation time sample resolution (seconds)
sim_smp = sim_bnd(2)/sim_res; % Initialize simulation time sample quanta (samples)
time = linspace(sim_res, sim_bnd(2), sim_smp); % Initialize sim time vector (seconds)

%% 2. MODEL PARAMETERS (PICK YOUR PREFERRED NEURON MODEL)

n = 16; % Set Izhikevich model parameters

p = [ 00.020,  00.200, -65.000,  06.000,  14.000;... % 01. Tonic spiking
      00.020,  00.250, -65.000,  06.000,  00.500;... % 02. Phasic spiking
      00.020,  00.200, -50.000,  02.000,  15.000;... % 03. Tonic bursting
      00.020,  00.250, -55.000,  00.050,  00.600;... % 04. Phasic bursting
      00.020,  00.200, -55.000,  04.000,  10.000;... % 05. Mixed-mode
      00.010,  00.200, -65.000,  08.000,  30.000;... % 06. Spike frequency adaptation
      00.020, -00.100, -55.000,  06.000,  00.000;... % 07. Class 1
      00.200,  00.260, -65.000,  00.000,  00.000;... % 08. Class 2
      00.020,  00.200, -65.000,  06.000,  07.000;... % 09. Spike latency
      00.050,  00.260, -60.000,  00.000,  00.000;... % 10. Subthreshold oscillations
      00.100,  00.260, -60.000, -01.000,  00.000;... % 11. Resonator
      00.020, -00.100, -55.000,  06.000,  00.000;... % 12. Integrator
      00.030,  00.250, -60.000,  04.000,  00.000;... % 13. Rebound spike
      00.030,  00.250, -52.000,  00.000,  00.000;... % 14. Rebound burst
      00.030,  00.250, -60.000,  04.000,  00.000;... % 15. Threshold variability
      01.000,  01.500, -60.000,  00.000, -65.000;... % 16. Bistability
      01.000,  00.200, -60.000, -21.000,  00.000;... % 17. DAP
      00.020,  01.000, -55.000,  04.000,  00.000;... % 18. Accomodation
     -00.020, -01.000, -60.000,  08.000,  80.000;... % 19. Inhibition-induced spiking
     -00.026, -01.000, -45.000,  00.000,  80.000];   % 20. Inhibition-induced bursting

%% 3. MEMBRANE TIME PARAMETERS

V_time = zeros(size(time)); % Allocate membrane potential time vector memory
U_time = V_time;            % Allocate membrane recovery variable time vector memory
V = -70;                    % Initial membrane potential
U = p(n, 2)*V;              % Initial membrane recovery
spike = V_time;             % Allocate spike vector memory
spike_time = V_time;        % Allocate spike time vector memory
spike_freq =  V_time;       % Allocate spike time vector memory

%% 4. STIMULUS

%  A. STIMULUS TIME SERIES (ALTER stim_bnd TO ALTER INJECTION CURRENT START/END)
stim_bnd = [0, 1]; % Initialize forcing signal generation time (seconds)
stim_int = diff(stim_bnd); % Initialize forcing signal time interval (seconds)
stim_smp = stim_int/sim_res; % Initialize forcing signal time sample quanta (samples)
stim = linspace(sim_res, stim_int, stim_smp); % Initialize forcing signal time vector values

%  B. STIMULUS AMPLITUDE (PLAY AROUND WITH amp AND amp_damp)
amp = [1, 1]; % Initialize carrier signal amplitude (mA)
amp_damp = 0; % Initialize terminal carrier signal dampening (mA)
% (BELOW) Initialize carrier signal amplitude vector values (mA)
amp = (linspace(amp(1), amp(2), stim_smp) + p(n, 5)).*exp(amp_damp*stim);

%  C. STIMULUS FREQUENCY (PLAY AROUND WITH freq_bnd AND freq_damp)
freq_bnd = [1, 1]; % Initialize carrier signal frequency (Hz)
freq_damp = 0; % Initialize terminal carrier signal dampening (Hz)
% (BELOW) Initialize carrier signal frequency vector values (Hz)
freq = (2*pi*linspace(freq_bnd(1), freq_bnd(2), stim_smp)).*exp(freq_damp*stim);;

%  D. STIMULUS PHASE (PLAY AROUND WITH phase_offset)
phase_offset = (2*pi/360)*0; % Convert initial forcing signal phase offset (degrees)
% Initialize carrier signal phase vector values (degrees)
phase = cumsum(freq)/sim_res; % Initialize carrier signal phase vector values (degrees)

%  E. STIMULUS NOISE (PLAY AROUND WITH amp_noise AND freq_noise)
% (BELOW) Initialize carrier signal amplitude noise vector values (percent)
amp_noise = amp.*randn(1, stim_smp)*0.01;
% (BELOW) Initialize carrier signal frequency noise vector values (percent)
freq_noise = freq.*randn(1, stim_smp)*0.01;

%  F. STIMULUS WAVEFORM GENERATION
stim = (amp + amp_noise).*sin((freq + freq_noise).*stim + phase_offset);
stim = [zeros(1, stim_bnd(1)/sim_res), stim, zeros(1, (sim_bnd(2) - stim_bnd(2))/sim_res)];

%% 5. SIMULATION

%  A. EXECUTION
for t = 1 : length(time) - 1
    
    V = V + 0.25*(0.04*V^2 + 5*V + 140 - U + stim(t));
    U = U + 0.25*p(n, 1)*(p(n, 2)*V - U);
    
    if V > 30                   % If there is a spike
       V_time(t + 1) = 30;      % Generate peak
       V = p(n, 3);             % Reset membrane voltage to c
       U = U + p(n, 4);         % Membrane recovery variable increment
       spike(t) = 1;            % Record spike
       spike_time(t) = time(t); % Record spike time
       spike_stim(t) = stim(t);
    else
        V_time(t + 1) = V;
    end;
    U_time(t + 1) = U;
end;

%  B. ANALYSIS
% (BELOW) Calculate spike train frequency
spike_freq = [0 sim_smp.*diff(spike)];
spike_freq = cumsum(spike_freq)/sim_res;
% (BELOW) Convolve spike train using gaussian kernel
sigma = 1;
x = [-0.1: 0.01 : 0.1];
k = exp(-(x/sigma).^2/2)/(sigma*sqrt(2*pi));
spike_freq = conv(spike_freq, k);
spike_freq = spike_freq(floor(length(k)/2) : end - ceil(length(k)/2));

%% OUTPUT

% (BELOW) Plot simulation variables
figure();
subplot(5, 1, 1); plot(time, stim);
title('Stimulus'); xlabel('Time (s)'); ylabel('mA');
subplot(5, 1, 2); plot(time, V_time);
title('Membrane Potential'); xlabel('Time (s)'); ylabel('mV');
subplot(5, 1, 3); plot(time, U_time);
title('Membrane Recovery Variable'); xlabel('Time (s)'); ylabel('mV');
subplot(5, 1, 4); plot(time, spike);
title('Spike Times'); xlabel('Time (s)'); ylabel('0/1');
subplot(5, 1, 5); plot(time, spike_freq, 'r');
title('Spike Frequency'); xlabel('Time (s)'); ylabel('Hz');