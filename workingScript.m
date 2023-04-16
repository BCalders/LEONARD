%% Doppler Localization

%% init
clear
clc
close all
format compact
format long

%% Setup
% global SAT RECIEVER C

simTime = 5;
startTime = datetime("5-july-2022 13:17");
stopTime = startTime + minutes(simTime);
sampleTime = 60;        % has to be 60 to be compliant with function

C = physconst("Lightspeed");

sc = satelliteScenario(startTime, stopTime, sampleTime);
% sc.AutoShow = false;

gs = groundStation(sc, 51.17800903509613, 4.418814450257098, 'Name', "CGB - Receiver");
gsEcefPos = lla2ecef([gs.Latitude, gs.Longitude, gs.Altitude])';

SAT.all = satellite(sc, "tle/iridium.tle");     % Iridium satellites used as a testing satellite set with global coverage
numSats = length(SAT.all);
SAT.femit = 1610e6;        % Avg emitted frequency in Hz used by Iridium
satAcc = repmat([0.9, 0.9, 0.9]', [1, numSats]);
initState = [4.6e+06, 1e+06, 4.2e+06, 0, 0];  % init pos in rome
%[4e+06, 3e+05, 5e+06, 0, 0]; init pos in the netherlands

%% Calculate all dopplerShifts for all timepoints

[dopV, fobs, r, vel] = dopShift(startTime, stopTime, SAT.all, gs, SAT.femit);

%% Calculation
% Define satellite position
figure
geoscatter(gs.Latitude, gs.Longitude, 'filled', 'MarkerFaceColor', 'r')
title("Initial Position Estimations")
hold on
llaState = ecef2lla(initState(1, 1:3));
geoscatter(llaState(1), llaState(2), 'filled', 'MarkerFaceColor','c')

llaState = [0,0,0]; % intialization

% define access for all timepoints
ac = access(SAT.all, gs);
acStatus = accessStatus(ac);

estimatedState = initState;

for currTime = 1:simTime+1
    focussedSat = 1;  % know which satellite of the satellites in view is being focussed on
    for currSat = 1:numSats
        if acStatus(currSat, currTime) == 1         % only calculate if satellite is in view

            % determine states
            [satPos,satVel] = states(SAT.all(currSat), startTime + minutes(currTime-1), "Coordinateframe", "ecef");
            satPos = squeeze(satPos);
            satVel = squeeze(satVel);

            % calculate vel and set previous vel
            if focussedSat > 1
                gs2satVelPrev = gs2satVel;
            else
                gs2satVelPrev = 1e+04;   % arbitrarily chosen value -> tbd
            end
            gs2satVel = calcRelVel(satPos, satVel, gsEcefPos);
            
            % calculate acceleration
            satAcc = gs2satVel - gs2satVelPrev;
            
            % calculate doppler shift
            fobs = speed2Dop(SAT.femit, gs2satVel);       % observed Doppler shift (Hz)
            
            % start calculation
            rangeVect = satPos - estimatedState(1:3)';
            rho = vecnorm(rangeVect);
            unitVector = rangeVect ./ rho;
            
            rhoDot = sum(satVel .* unitVector);
            rhoDotDot = sum(satAcc .* unitVector);
            
            eDot = (1.0 ./ rho) .* (satVel - unitVector .* rhoDot);

            if focussedSat > 1     % if it is the first timepoint, there is no H yet for this timepoint
                H = [H ; eDot', 1, -rhoDotDot'];
            else 
                H = [eDot', 1, -rhoDotDot'];
            end 
            
            relVel = dot(satVel, unitVector);                   % relative velocity (m/s)
            D_predicted = speed2Dop(SAT.femit, relVel);       % predicted Doppler shift (Hz));
            deltaDoppler = fobs(1) - D_predicted;           % deltaDoppler (Hz)
            deltaD = dop2Speed(SAT.femit, deltaDoppler);      % deltaD (m/s)
            
            if focussedSat > 1     % if it is the first satellite, there is no deltaDMatrix yet for this timepoint
                deltaDMatrix = [deltaDMatrix; deltaD];
            else
                deltaDMatrix = deltaD;
            end
            focussedSat = focussedSat + 1; % indicate that next satellite is being focussed on
        end
    end
    deltaX = H \ deltaDMatrix;                                % use of backslash operator to invert H

    estimatedState = estimatedState + deltaX';
    llaState(currTime, :) = ecef2lla(estimatedState(1:3));
    % drawnow;
    disp("Current position: x: " + estimatedState(1) + " y: " + estimatedState(2) + " z: " + estimatedState(3))
end
geoscatter(llaState(:, 1), llaState(:, 2), 'b')
figure;
geoscatter(gs.Latitude, gs.Longitude, 'filled', 'MarkerFaceColor', 'r')
hold on
geoscatter(llaState(:, 1), llaState(:, 2), 'b')
title("Position Estimations zoom")
geolimits([51.170833 51.1875], [4.4 4.433333])

% play(sc)  % run using F9 to show satelliteScenario

%% functions

function [fObserved] = speed2Dop(fEmit, relativeVelocity)
    vReceiver = 0;
    vSource = relativeVelocity;
    c = physconst("Lightspeed");
    fObserved = fEmit * ((c + vReceiver) ./ (c + vSource));
end

function [relativeVelocity] = dop2Speed(fEmit, fObserved)
    c = physconst("Lightspeed");
    relativeVelocity = c * (fObserved - fEmit) / fEmit;
end

function [relativeVelocity] = calcRelVel(satPos, satVel, recPos)
    unitVector = (satPos - recPos) ./ vecnorm(satPos - recPos);

    relativeVelocity = dot(satVel, unitVector);
end