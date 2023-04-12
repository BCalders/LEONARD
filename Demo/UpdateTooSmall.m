%% Doppler Localization

%% init
clear
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
SAT.fcarrier = 1610e6;        % Avg emitted frequency in Hz used by Iridium
% fo = 0;
[satPos,satVel] = states(SAT.all, startTime);
satAcc = repmat([0.9, 0.9, 0.9]', [1, numSats]);
initState = repmat([4e+06, 3e+05, 5e+06, 0, 0], [numSats,1]);

[dopV, fo, r, vel] = dopShift(startTime, stopTime, SAT.all, RECIEVER.gs, SAT.fcarrier);

%% Calculation
% Define satellite position
figure
geoscatter(RECIEVER.gs.Latitude, RECIEVER.gs.Longitude, 'filled', 'MarkerFaceColor', 'r')
title("Initial Position Estimations")
hold on

for currTime = 1:100 %simTime+1
temp = 1; % for debug purposes, has to be changed to currtime in the code again

[satPos,satVel] = states(SAT.all, startTime); % + seconds(currTime-1));
satPos = squeeze(satPos);
satVel = squeeze(satVel);

    for currSat = 1:numSats
        
        e = satPos(:, currSat) - initState(currSat,1:3);      
        rho = vecnorm(e')';                             
        e = e ./ rho;                                   
        
        rhoDot = sum((satVel(:, currSat) .* e),2);      
        rhoDotDot = sum((satAcc(:, currSat) .* e),2);  
        
        eDot = rhoDot./rho;                             
        
        H = [eDot, ones(3, 1), -rhoDotDot];            
        
        D_predicted = -sum((e .* satVel(:, currSat)),2);
        D_measured = fo(:, currSat) * C / SAT.fcarrier;                      % 
        deltaD = D_measured' - D_predicted';                                 % deltaD (m/s)    
        
        deltaX = H / deltaD;                            
        
        deltaRho = H*deltaX;                            
        
        initState(currSat, 1:3) = initState(currSat, 1:3) + deltaX'; 
        llaState = ecef2lla(initState(currSat, 1:3));
        geoscatter(llaState(1), llaState(2), 'b')
%         hold on
        disp("Current position: " + initState(currSat, 1:3))
    end

end

% geoscatter(RECIEVER.gs.Latitude, RECIEVER.gs.Longitude, 'filled', 'MarkerFaceColor', 'r')
% title("Initial Position Estimations")


% verjlaren van de variabelen
% dopshift herwerken
% focus op de werking niet efficientie
% assume locatie satellieten gekend
% 1 maart pitches -> picht moet zijn van level zodat de masterstudenten kunnen
% volgen