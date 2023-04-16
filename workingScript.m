%% Doppler Localization

%% init
clear
clc
close all
format compact
format long

%% Setup
% global SAT RECIEVER C

simTime = 10;
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
initState = [4e+06, 3e+05, 5e+06, 0, 0]; %repmat([4e+06, 3e+05, 5e+06, 0, 0], [numSats,1]);

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

% define access for all timepoints
ac = access(SAT.all, gs);
acStatus = accessStatus(ac);

for currTime = 1:simTime+1



    % filter out satellites that are not in view
    allSats = SAT.all(acStatus(:, currTime));

    % determine states of all satellites in view
    [satPos,satVel] = states(allSats, startTime + minutes(currTime-1), "Coordinateframe", "ecef");
    numSats = length(satPos);
    satPos = squeeze(satPos);
    satVel = squeeze(satVel);

    % calculate acceleration of satellites in view
    % if currTime > 1     % if it is the first timepoint, there is no previous timepoint to calculate acceleration
    %     satAcc = satVel - satVelPrev;
    % else
    %     satAcc = repmat([0.9, 0.9, 0.9]', [1, numSats]);
    % end

    % calculate acceleration of satellites in view
    if currTime > 1
        gs2satVelPrev = gs2satVel;
    else
        gs2satVelPrev = zeros(1, numSats);
    end
    gs2satVel = calcRelVel(satPos, satVel, gsEcefPos);
    
    satAcc = gs2satVel - gs2satVelPrev;
    
    % calculate doppler shift for all satellites in view
    fobs = speed2Dop(SAT.femit, gs2satVel);       % observed Doppler shift (Hz)

    rangeVect = satPos - initState(1:3)';            % e(XYZvalue, satNr)
    rho = vecnorm(rangeVect);                             
    e = rangeVect ./ rho;                                   

    rhoDot = sum(satVel .* e);      
    rhoDotDot = sum(satAcc .* e);  

    eDot = (1.0 ./ rho) .* (satVel - e .* rhoDot);
    % eDot = gradient(e);
    % eDot = (rhoDot - e * dot(rho, rhoDot, 1))/rho^2;

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
% play(sc)  % run using F9 to show satelliteScenario

%% functions

function [fObserved] = speed2Dop(fEmit, relativeVelocity)
    vReceiver = 0;
    vSource = relativeVelocity;
    c = physconst("Lightspeed");
    fObserved = fEmit * ((c + vReceiver) ./ (c + vSource));
    % fObserved = fEmit * (1 + relativeVelocity / C); % TODO c+ vr/ c+vs juiste implementeren
end

function [relativeVelocity] = dop2Speed(fEmit, fObserved)
    c = physconst("Lightspeed");
    relativeVelocity = c * (fObserved - fEmit) / fEmit;
end

function [relativeVelocity] = calcRelVel(satPos, satVel, recPos)
    unitVector = (satPos - recPos) ./ vecnorm(satPos - recPos);

    relativeVelocity = dot(satVel, unitVector);
end


%% to be solved
    % make sure the acc is calculated using the right satellites 
    % 