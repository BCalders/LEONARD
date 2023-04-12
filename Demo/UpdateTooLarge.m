%% Doppler Localization

%% init
clear
clc
close all
format compact
format long

%% Setup
% global SAT RECIEVER C

simTime = 1;
startTime = datetime("5-july-2022 13:17");
stopTime = startTime + minutes(simTime);
sampleTime = 60;        % has to be 60 to be compliant with function

C = physconst("Lightspeed");

sc = satelliteScenario(startTime, stopTime, sampleTime);
sc.AutoShow = false;

RECIEVER.gs = groundStation(sc, 51.17800903509613, 4.418814450257098, 'Name', "CGB - Receiver");
RECIEVER.pos = lla2ecef([RECIEVER.gs.Latitude, RECIEVER.gs.Longitude, RECIEVER.gs.Altitude]);

SAT.all = satellite(sc, "../tle/iridiumFilter.tle");     % Iridium satellites used as a testing satellite set with global coverage
numSats = length(SAT.all);
SAT.femit = 1610e6;        % Avg emitted frequency in Hz used by Iridium
satAcc = repmat([0.9, 0.9, 0.9]', [1, numSats]);
initState = [4e+06, 3e+05, 5e+06, 0, 0]; %repmat([4e+06, 3e+05, 5e+06, 0, 0], [numSats,1]);

[dopV, fobs, r, vel] = dopShift(startTime, stopTime, SAT.all, RECIEVER.gs, SAT.femit);

%% Calculation
% Define satellite position
figure
geoscatter(RECIEVER.gs.Latitude, RECIEVER.gs.Longitude, 'filled', 'MarkerFaceColor', 'r')
title("Initial Position Estimations")
hold on
llaState = ecef2lla(initState(1, 1:3));
geoscatter(llaState(1), llaState(2), 'filled', 'MarkerFaceColor','c')

for currTime = 1:1 %simTime+1
temp = 1; % for debug purposes, has to be changed to currtime in the code again

[satPos,satVel] = states(SAT.all, startTime); % + seconds(currTime-1));
satPos = squeeze(satPos);
satVel = squeeze(satVel);

% TODO implement as a for loop to have less errors

    e = satPos - initState(1:3)';            % e(XYZvalue, satNr)
    rho = vecnorm(e);                             
    e = e ./ rho;                                   
    
    rhoDot = sum(satVel .* e);      
    rhoDotDot = sum(satAcc .* e);  
    
    eDot =  (1.0 ./ rho) .* (satVel - e .* rhoDot);
    
    H = [eDot', ones(6, 1), -rhoDotDot'];            
    
    relVel = dot(satVel, e);                   % relative velocity (m/s)
    D_predicted = speed2Dop(SAT.femit, relVel);       % predicted Doppler shift (Hz));
    deltaDoppler = fobs(1) - D_predicted;           % deltaDoppler (Hz)
    deltaD = dop2Speed(SAT.femit, deltaDoppler);      % deltaD (m/s)
    
    deltaX = H \ deltaD';                                % use of backslash operator to invert H
                             
    
    initState = initState + deltaX'; 
    llaState = ecef2lla(initState(1:3));
    geoscatter(llaState(1), llaState(2), 'b')
    disp("Current position: x: " + initState(1) + " y: " + initState(2) + " z: " + initState(3))
end

%% functions

function [fObserved] = speed2Dop(fEmit, relativeVelocity)
    C = physconst("Lightspeed");
    fObserved = fEmit * (1 + relativeVelocity / C); % TODO c+ vr/ c+vs juiste implementeren
end

function [relativeVelocity] = dop2Speed(fEmit, fObserved)
    C = physconst("Lightspeed");
    relativeVelocity = C * (fObserved - fEmit) / fEmit;
end