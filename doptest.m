% Coarse-time Doppler navigation equations for a static receiver based on Snapshot positioning without initial information
% by Fernández-Hernández I and Borre K 
% implemented by Bas Calders 2023 
% version 1.1

clear 
clc
close all

% Define constants
c = physconst("Lightspeed"); % Speed of light in m/s
f_L1 = 1575.42e6; % L1 frequency in Hz
lambda_L1 = c/f_L1; % L1 wavelength in m
f_L2 = 1227.60e6; % L2 frequency in Hz
lambda_L2 = c/f_L2; % L2 wavelength in m
f_IF = 9.43e6; % Intermediate frequency in Hz
lambda_IF = c/f_IF; % Intermediate frequency wavelength in m
t_s = 1/50; % Sampling period in seconds
N = 4; % Number of satellites

% Define satellite positions
satPos = [15600e3, 7540e3, 20140e3;
    18760e3, 2750e3, 18410e3;          
    19170e3, 14690e3, 13180e3;          
    26560e3, 1230e3, 11270e3];


% Generate simulated signal
t = (0: t_s : 1-t_s); % Time vector
phi = 2*pi*f_L1*t; % Carrier phase
c_a = (1+cos(phi))/2; % C/A code
c_p = sin(2*phi); % P code
codeDelay = [0, 1, 2, 3]; % Code delays for each satellite
carrierFreq = [f_L1, f_L2, f_L1, f_L2]; % Carrier frequencies for each satellite
carrierPhase = [0, 0, 0, 0]; % Carrier phases for each satellite
for i = 1:N
    signal = signal + (c_a.*cos(2*pi*carrierFreq(i)*t + carrierPhase(i) + codeDelay(i)*lambda_L1/lambda_IF) + ...
        c_p.*cos(2*pi*carrierFreq(i+1)*t + carrierPhase(i+1) + ...
        codeDelay(i)*lambda_L2/lambda_IF)).*exp(-1i*2*pi*sqrt((satPos(i,1))^2 + ...
        (satPos(i,2))^2 + (satPos(i,3))^2)/lambda_L1);
end 

% Coarse-time processing
cpr = zeros(1,N); % Initialize code phase range
for i = 1:N
    [~,ind] = max(abs(xcorr(signal, c_a.*cos(2*pi*carrierFreq(i)*t + carrierPhase(i)), 'normalized'))); % Find peak of C/A code correlation
    cpr(i) = mod(ind - 1 - codeDelay(i)*t_s*f_L1, length(c_a)); % Calculate code phase range
end

% Calculate pseudoranges
pr = zeros(1,N); % Initialize pseudoranges
for i = 1:N
    pr(i) = cpr(i)*lambda_L1 + satPos(i,1); % Calculate pseudorange
end

% Solve for receiver position
A = zeros(N-1, 3); % Initialize A matrix
b = zeros(N-1, 1); % Initialize b vector
for i = 1:N-1
    A(i,:) = (satPos(i+1,:) - satPos(1,:))/norm(satPos(i+1,:) - satPos(1,:)) - (satPos(i+1,:) - satPos(i,:))/norm(satPos(i+1,:) - satPos(i,:)); % Populate A matrix
    b(i) = pr(i+1) - pr(1); % Populate b vector
end
x = A\b; % Solve for receiver position

% Display receiver position
fprintf('Receiver position: (%.4f, %.4f, %.4f) meters\n', x(1), x(2), x(3));

