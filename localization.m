%% Doppler Localization

%% init
clear
% close all
format compact
format long

%% Setup
global SAT RECIEVER C

simTime = 1;
startTime = datetime("5-july-2022 13:17");
stopTime = startTime + minutes(simTime);
sampleTime = 60;        % has to be 60 to be compliant with function

C = physconst("Lightspeed");

sc = satelliteScenario(startTime, stopTime, sampleTime);
sc.AutoShow = false;

% receiver = 'CGB';       % Adding an arbitrary ground based reciever for testing purposes

% switch(receiver)
%     case 'CGB'
%         RECIEVER.gs = groundStation(sc, 51.17800903509613, 4.418814450257098, 'Name', "CGB - Receiver");
%     case 'UC'
%         RECIEVER.gs = groundStation(sc, 37.871946249596625, -122.25853766615649, 'Name', "UC - Receiver");
%     case 'USYD'
%         RECIEVER.gs = groundStation(sc, -33.88857476158162, 151.1873333064266, 'Name', "USYD - Receiver");
%     case 'UBA'
%         RECIEVER.gs = groundStation(sc, -34.59978022088964, -58.373369858300805, 'Name', "UBA - Receiver");
%     otherwise
%         RECIEVER.gs = groundStation(sc, 51.17800903509613, 4.418814450257098, 'Name', "CGB - Receiver");
% end

RECIEVER.gs = groundStation(sc, 51.17800903509613, 4.418814450257098, 'Name', "CGB - Receiver");
RECIEVER.pos = lla2ecef([RECIEVER.gs.Latitude, RECIEVER.gs.Longitude, RECIEVER.gs.Altitude]);


SAT.all = satellite(sc, "tle/iridiumFilter.tle");     % Iridium satellites used as a testing satellite set with global coverage
numSats = length(SAT.all);
SAT.fcarrier = 1610e6;        % Avg emitted frequency in Hz used by Iridium
% fo = 0;
disp('Setup complete.')

%% Calculating doppler
disp('Calculating Doppler shift...')
[dopv, dopFreq, r, vel] = dopShift(startTime, stopTime, SAT.all, RECIEVER.gs, SAT.fcarrier);

%% Calculation
% Define satellite positions
satPos = [15600e3, 7540e3, 20140e3;
    18760e3, 2750e3, 18410e3;          
    19170e3, 14690e3, 13180e3;          
    26560e3, 1230e3, 11270e3];

[satPos, satVel, satAcc] = get_sat_pos_vel_acc(SAT.all);

e = satPos - initState(1:3);
rho = vecnorm(e')';
e = e ./ rho;

rhoDot = sum((satVel .* e)')';
rhoDotDot = sum((satAcc .* e)')';

H = [eDot, ones(length(acqSats), 1), -rhoDotDot];

D_predicted = -sum((e .* satVel)')' + initState(4);
D_measured = dopv' .* C ./ SAT.fcarrier;
deltaD = D_measured' - D_predicted;

deltaX = H \ deltaD;

deltaRho = H*deltaX; % assuming no noise deltaRho = update of the a priori state


%% get sat pos vel and acc

function [satPos, satVel, satAcc] = get_sat_pos_vel_acc(sat)
    [satPos,satVel] = states(sat);

end

% verjlaren van de variabelen
% dopshift herwerken
% focus op de werking niet efficientie
% assume locatie satellieten gekend
% 1 maart pitches -> picht moet zijn van level zodat de masterstudenten kunnen
% volgen