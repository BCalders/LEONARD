%% Doppler Localization

%% Notes - Information

% After running 3 different colors of dots are visible on the map. the Cyan
% dot is the initial estimation. At this point it is set to Rome for added
% difficulty. The blue dots are the estimated positions after each iteration 
% they converge to the actual position of the ground station. The red dot
% is the actual position of the ground station.
% As a proof of concept the Starlink satellites are also avaliable in the
% repo. But please note this takes a long time to run.

%% init
clear
clc
close all
format compact
format long

%% Setup

disp("Setting up...")

simTime = 10;
% startTime = datetime("5-july-2022 13:17");
% startTime = datetime("11-september-2022 17:53");
% startTime = datetime("7-march-2023 04:22");
startTime = datetime("16-april-2023 16:05:33");
stopTime = startTime + minutes(simTime);
sampleTime = 10;
samplesPerMinute = sampleTime/60;

C = physconst("Lightspeed");

sc = satelliteScenario(startTime, stopTime, sampleTime);

gs = groundStation(sc, 51.17800903509613, 4.418814450257098, 'Name', "CGB - Receiver");
gsEcefPos = lla2ecef([gs.Latitude, gs.Longitude, gs.Altitude])';

SAT.all = satellite(sc, "tle/iridium.tle");     % Iridium satellites used as a testing satellite set with global coverage
numSats = length(SAT.all);
SAT.femit = 1610e6;        % Avg emitted frequency in Hz used by Iridium
satAcc = repmat([0.9, 0.9, 0.9]', [1, numSats]);
initState = [4.6e+06, 1e+06, 4.2e+06, 0, 0];  % init pos in Rome for added difficulty
%[4e+06, 3e+05, 5e+06, 0, 0]; init pos in the netherlands

breakLoop = false;
loopBroken = false;

disp("Setup complete")

%% Calculation

disp("Starting calculation...")

% Define satellite position
figure
geoscatter(gs.Latitude, gs.Longitude, 'filled', 'MarkerFaceColor', 'r')
title("Initial Position Estimations")
hold on
llaStates = ecef2lla(initState(1, 1:3));
geoscatter(llaStates(1), llaStates(2), 'filled', 'MarkerFaceColor','c')

% define access for all timepoints
ac = access(SAT.all, gs);
acStatus = accessStatus(ac);

estimatedState = initState;

for currTime = 1:simTime
    focussedSat = 1;                                % know which satellite in view is being focussed on
    for currSat = 1:numSats
        if acStatus(currSat, currTime) == 1         % only calculate if satellite is in view

            % determine satellite position and velocity
            [satPos,satVel] = states(SAT.all(currSat), startTime + minutes(currTime), "Coordinateframe", "ecef"); % using currtime and not currtime - 1 to ensure a previous time is present in the SC 
            satPos = squeeze(satPos);
            satVel = squeeze(satVel);

            % setting previous velocity
            if focussedSat > 1
                gs2satVelPrev = gs2satVel;
            else
                % [satPos,satVel] = states(SAT.all(currSat), startTime + minutes(currTime-2), "Coordinateframe", "ecef");
                % satPos = squeeze(satPos);
                % satVel = squeeze(satVel);
                % gs2satVelPrev = calcRelVel(satPos, satVel, gsEcefPos); 
                gs2satVelPrev = 1e+04;   % arbitrarily chosen value -> tbd
            end

            % calculate relative velocity
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

%             if deltaD == -C                             % break loops if accuracy can not become any better
%                 disp("maximum achievable accuracy is highly likely") % TODO -> move this to time loop -> deltaDMatrix should be all C for loopbreak
%                 if breakLoop
%                     loopBroken = true;
%                 end
%                 breakLoop = true;
%             end

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
    results(currTime,:) = estimatedState;
    llaStates(currTime, :) = ecef2lla(estimatedState(1:3));
    % drawnow;
    disp(currTime + ": Current position: x: " + estimatedState(1) + " y: " + estimatedState(2) + " z: " + estimatedState(3))
    
    if loopBroken               % break loop when accuracy is perfect -> save calculations
        disp("maximum achievable accuracy reached")
        break;
    end
end
disp("Calculation complete")

disp("plotting...")

geoscatter(llaStates(:, 1), llaStates(:, 2), 'b')
figure;
geoscatter(gs.Latitude, gs.Longitude, 'filled', 'MarkerFaceColor', 'r')
hold on
geoscatter(llaStates(:, 1), llaStates(:, 2), 'xb')
title("Position Estimations zoom")
geolimits([51.170833 51.1875], [4.4 4.433333])
% geobasemap('streets')

%% plot results


for i = 1:size(results, 1)
    distError(i) = vecnorm(gsEcefPos - results(i, 1:3)');
end
figure
plot(distError)

disp("I'm done!")

%% functions

% play(sc)  % run using F9 to show satelliteScenario

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

%% Notes - To Do

% - First Acc value is now just hardcoded, should be derived correctly
% - Make an error plot at different ti
% mes - deviations/number of satellites used 
% - Attempt this implementation for moving User
% - First drafts for moving implementation -> to see differences with
% performance of this algorithm
% x Sometimes the first iteration is at coordinate 0 0 -> to be
% investigated => could be result of first acc estimation
% x completely phased out dopShift function
% - update Readme.md file to better image the current state of the
% algorithm
% x less time between measurements
% x resultaat -> accuracy vs tijd
% - nieuwe UE locatie -> bestemming afrika ;)
% - better way to break loop


%% Notes - Solutions

% -
% -
% -
% -
% - foutje bij het schrijven van resultaten -> timings kwamen niet correct
% overeen met llaState waarden 
% - ok
% - 
% - hard to fix because I work in minutes not seconds
% - first implementation done
% -
% -