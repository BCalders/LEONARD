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
startTime = datetime("7-march-2023 04:22:00");
% startTime = datetime("16-april-2023 16:40:30");
stopTime = startTime + seconds(simTime);
sampleTime = 1;

c = physconst("Lightspeed");

sc = satelliteScenario(startTime, stopTime, sampleTime);
% add even more locations (every continent)
% gs = groundStation(sc, 51.17800903509613, 4.418814450257098, 'Name', "CGB - Receiver");
% gs = groundStation(sc, 0.5108574230657834, 33.13331679803374, 'Name', "Uganda - Reciever");
gs = groundStation(sc, -12.051334463667322, -77.012949622246, 'Name', "Peru - Reciever");
gsEcefPos = lla2ecef([gs.Latitude, gs.Longitude, gs.Altitude])';

% globalStar
SAT.all = satellite(sc, "tle/globalstar.tle");     % GlobalStar satellites used as a testing satellite set with global coverage 
numSats = length(SAT.all);
SAT.femit = 16118.725e6;        % Avg emitted frequency in Hz used by GlobalStar in L-band for mobile ue
% initState = [4.6e+06, 1e+06, 4.2e+06, 2.99, 100];  % init pos in Rome for added difficulty
initState = [4.12e+06, -4.55e+06, -1.72e+06, 0, 0];  % init pos in Brasilia for added difficulty
% initState = [6.13e+06, 1.68e+06, -4.9e+05, 0, 0];  % init pos in Kinshasa for added difficulty
%[4e+06, 3e+05, 5e+06, 0, 0]; init pos in the netherlands

disp("Setup complete")

figure('Name', "Map Plot")
geoscatter(gs.Latitude, gs.Longitude, 'filled', 'MarkerFaceColor', 'r')
title("Initial Position Estimations")
hold on
llaStates = ecef2lla(initState(1, 1:3));
geoscatter(llaStates(1), llaStates(2), 'filled', 'MarkerFaceColor','c')

%% Calculation

disp("Starting calculation...")

% define access for all timepoints
ac = access(SAT.all, gs);
acStatus = accessStatus(ac);

distError = zeros(simTime, 1);


estimatedState = initState;

for currTime = 1:simTime
    focussedSat = 0;                                % know which satellite in view is being focussed on
    satsInView = nnz(acStatus(: ,currTime));

    if satsInView < 4
        disp("Not Enough satellites to perform localization")
        break;
    end

    deltaD = zeros(satsInView, 1);
    H = zeros(satsInView, 5);

    for currSat = 1:numSats
        if acStatus(currSat, currTime) == 1         % only calculate if satellite is in view

            focussedSat = focussedSat + 1;          % indicate that next satellite is being focussed on

            % determine satellite position and velocity
            [satPos,satVel] = states(SAT.all(currSat), startTime + seconds(currTime), "Coordinateframe", "ecef");  % using currtime and not currtime - 1 to ensure a previous time is present in the SC 

            % setting previous velocity
            [satPosPrev,satVelPrev] = states(SAT.all(currSat), startTime + seconds(currTime - 1), "Coordinateframe", "ecef");
            sat2gsVelPrev = calcRelVel(satPosPrev, satVelPrev, gsEcefPos); 

            % calculate relative velocity
            sat2gsVel = calcRelVel(satPos, satVel, gsEcefPos);  % relative velocity using TLE info(m/s)
            
            % calculate acceleration
            satAcc = (satVelPrev - satVel)./sampleTime;
            
            % calculate doppler shift
            fobs = speed2Dop(SAT.femit, sat2gsVel);       % observed Doppler shift (Hz)
            
            % start calculation
            rangeVect = satPos - estimatedState(1:3)';
            rho = vecnorm(rangeVect);
            unitVector = rangeVect / rho;
            
            rhoDot = dot(satVel, unitVector);               % Range Rate (m/s)
            rhoDotDot = dot(satAcc, unitVector);            % Range Acceleration (m/s^2)
            
            eDot = (1/rho) * (satVel - unitVector * rhoDot);    % Unit Vector Derivative

            H(focussedSat,:) = [eDot', 1, -rhoDotDot];          % Jacobian Matrix
            
            D_predicted = rhoDot - estimatedState(4);           % predicted D (m/s)
            deltaD(focussedSat) = sat2gsVel - D_predicted;      % deltaD (m/s)
        end
    end
    deltaX = H \ deltaD;                                % use of backslash operator to invert H
    
    estimatedState = estimatedState - deltaX';
    results(currTime,:) = estimatedState;
    llaStates(currTime, :) = ecef2lla(estimatedState(1:3));
    % drawnow;
    disp(currTime + ": Current Estimated position: x: " + estimatedState(1) + " y: " + estimatedState(2) + " z: " + estimatedState(3))
end
disp("Calculation complete")

disp("plotting...")

names = string(1:simTime);  
geoscatter(llaStates(:, 1), llaStates(:, 2), 'b')
text(llaStates(:,1),llaStates(:,2),names)
figure('Name', "Zoomed")
geoscatter(gs.Latitude, gs.Longitude, 'filled', 'MarkerFaceColor', 'r')
hold on
geoscatter(llaStates(:, 1), llaStates(:, 2), 'xb')
text(llaStates(:,1),llaStates(:,2),names)
title("Position Estimations zoom")
geolimits([gs.Latitude - 0.001,  gs.Latitude + 0.001], [gs.Longitude - 0.02 gs.Longitude + 0.02])
% geobasemap('streets')

%% plot results


for i = 1:simTime
    distError(i) = vecnorm(gsEcefPos - results(i, 1:3)');
end

figure('Name', "Error Plot")
title("Distance Errors from second iteration")
plot(distError(1:simTime))

disp("I'm done!")

%% functions

% play(sc)  % run using F9 to show satelliteScenario

function [fObserved] = speed2Dop(fEmit, relativeVelocity)
    vReceiver = 0;
    vSource = relativeVelocity;
    c = physconst("Lightspeed");
    fObserved = fEmit * ((c + vReceiver) / (c + vSource));
end

function [relativeVelocity] = dop2Speed(fEmit, fObserved)
    c = physconst("Lightspeed");
    relativeVelocity = c * (fObserved) / fEmit;
end

function [relativeVelocity] = calcRelVel(satPos, satVel, recPos)
    unitVector = (satPos - recPos) ./ vecnorm(satPos - recPos);

    relativeVelocity = dot(satVel, unitVector);
end

%% Notes

% results: accuracy vs num of sats dependecy
% everything in excel -> parameters and results

% input parameters -> kind of constellation, init
% location(randomized/different distaces) 

% 15/5
% moving ue
% beschrijven van alles, sats hoogte, freq en alles algos dopp -> dan resultaten
% die hebben ook grafieken alle getallen. discussie is waarom uitschrijven
% resutaten en discussie in puntjes
% intro sota en conclusie kan gelijkaardig op vorige
% 23/5 
% moving ue zeker klaar
% eerste uitgeschreven versie doorsturen