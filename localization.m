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
[dopv, fo, r, vel] = dopShift(startTime, stopTime, SAT.all, RECIEVER.gs, SAT.fcarrier);

%% show vis

groundTrack(SAT.all,"LeadTime",1200);
sv = satelliteScenarioViewer(sc, "CameraReferenceFrame","Inertial", ...
    "Name","DopplerVis", "ShowDetails",true, ...
    "Basemap","colorterrain", "PlaybackSpeedMultiplier",0);
play(sv)

%% Calculation

[satPos, satVel, satAcc] = get_sat_pos_vel_acc(tow, ephFilt);
    
% e(satNr, 1) = (satPos(1) - initState(1)) / r;
% e(satNr, 2) = (satPos(2) - initState(2)) / r;
% e(satNr, 3) = (satPos(3) - initState(3)) / r;
e = satPos - initState(1:3);
rho = vecnorm(e')'; % normword of each row
e = e ./ rho; % [e1_x/rho1, e1_y/rho1, e1_z/rho1; e2_x/rho1, ...]

rhoDot = sum((satVel .* e)')';
rhoDotDot = sum((satAcc .* e)')';

eDot = (1.0 ./ rho) .* (satVel - e .* rhoDot);

H = [eDot, ones(length(acqSats), 1), -rhoDotDot];

D_predicted = -sum((e .* satVel)')' + initState(4);
D_measured = doppler' .* settings.c ./ 1575.42e6;
deltaD = D_measured' - D_predicted;

deltaX = H \ deltaD;

initState = initState + deltaX';
