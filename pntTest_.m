%% First Doppler PNT test

%% init
clear
% close all
format compact
format long

% addpath ekfukftoolbox\ %'autoarrange figures'

%% Setup
global SAT RECIEVER C

simTime = 0;
startTime = datetime("5-july-2022 13:17");
stopTime = startTime + minutes(simTime);
sampleTime = 60;        % has to be 60 to be compliant with function

C = physconst("Lightspeed");

sc = satelliteScenario(startTime, stopTime, sampleTime);
sc.AutoShow = false;

receiver = 'CGB';       % Adding an arbitrary ground based reciever for testing purposes

switch(receiver)
    case 'UC'
        RECIEVER.gs = groundStation(sc, 37.871946249596625, -122.25853766615649, 'Name', "UC - Receiver");
    case 'USYD'
        RECIEVER.gs = groundStation(sc, -33.88857476158162, 151.1873333064266, 'Name', "USYD - Receiver");
    case 'UBA'
        RECIEVER.gs = groundStation(sc, -34.59978022088964, -58.373369858300805, 'Name', "UBA - Receiver");
    case 'CGB'
    otherwise
        RECIEVER.gs = groundStation(sc, 51.17800903509613, 4.418814450257098, 'Name', "CGB - Receiver");
end
RECIEVER.pos = lla2ecef([RECIEVER.gs.Latitude, RECIEVER.gs.Longitude, RECIEVER.gs.Altitude]);

SAT.all = satellite(sc, "tle/iridiumFilter.tle");     % Iridium satellites used as a testing satellite set with global coverage
numSats = length(SAT.all);
SAT.fcarrier = 1610e6;        % Avg emitted frequency in Hz used by Iridium
% fo = 0;
disp('Setup complete.')

%% Dopshift Calc
disp('Calculating Doppler shift...')
tic
[dopv, fo, ~, vel] = dopShift(startTime, stopTime, SAT.all, RECIEVER.gs, SAT.fcarrier);
time = toc;
disp(['Calculation took: ', num2str(time), 's'])

usableSats = ~isnan(fo(simTime+1,:));
[spos, ~] = states(SAT.all, startTime + minutes(simTime), 'CoordinateFrame', 'ecef');
SAT.pos = squeeze(spos);

disp('Calculating Initial location guesses')
calc = tic;
coords = calcLocation(fo, usableSats, vel);
time = toc(calc);
disp(['Calculation took: ', num2str(time), 's'])
hold off
geoscatter(coords(1, :), coords(2, :), 'b')
Error = sqrt(sum((coords - [RECIEVER.gs.Latitude, RECIEVER.gs.Longitude, RECIEVER.gs.Altitude]').^2, 1,'omitnan'));
plot(Error)
xlim([60 106])
ylim([0 1e7])
title("Position Estimation Error")
xlabel("Estimation number")
ylabel("Error (deg)")
%% filter
disp('Starting Filter...')
[results, error, avgResult] = particlefilt(fo, startTime, stopTime, coords, sc);

%% Plotting estimations

figure
disp('Plotting results.')
geoscatter(coords(1, :), coords(2, :), 'b')
hold on
% geoscatter(avgResult(1), avgResult(2), 'filled', 'MarkerFaceColor', 'm')
geoscatter(RECIEVER.gs.Latitude, RECIEVER.gs.Longitude, 'filled', 'MarkerFaceColor', 'r')
% geobasemap('streets')
title("Initial Position Estimations")

disp('Goodbye!')


function coords = calcLocation(fo, usableSats, vel)
global SAT C
        p = sym('p',[1 3]);

        allsats = find(usableSats == 1);
        N = length(allsats);
        binary = dec2bin(0:2^N-1)' - '0';
        options = logical(binary(:, sum(binary, 1) == 3))';

        reverseStr = '';
        time = 0;

        for it = 1:length(options)
            remainder = tic;
            sats = allsats(options(it, :));

            delta = (C * (fo(:, sats)))/(SAT.fcarrier);
%             eqn = + delta.^2 .* ((p(1) - SAT.pos(1, sats)).^2 + (p(2) - SAT.pos(2, sats)).^2 + (p(3) - SAT.pos(3, sats)).^2) ...
%                 - (vel(1, sats) .* (p(1) + SAT.pos(1, sats)) +vel(2, sats) .* (p(2) + SAT.pos(2, sats)) + vel(3, sats) .* (p(3) - SAT.pos(3, sats))).^2 == 0;
            eqn = + delta.^2 .* ((SAT.pos(1, sats)-p(1)).^2 + (SAT.pos(2, sats)-p(2)).^2 + (SAT.pos(3, sats)-p(3)).^2) ...
                - (vel(1, sats) .* (SAT.pos(1, sats)-p(1)) + vel(2, sats) .* (SAT.pos(2, sats)-p(2)) + vel(3, sats) .* (SAT.pos(3, sats)-p(3))).^2 == 0;
            solution = solve(eqn, p);
            
            p1 = double(solution.p1);
            p2 = double(solution.p2);
            p3 = double(solution.p3);

            lla = ecef2lla(real([p1, p2, p3]));
            if exist('coords', 'var') == 1
                coords = [coords, lla'];
            else
                coords = lla';
            end
            time2 = toc(remainder);
            time = time + time2;
            remaining = time/it * (length(options) - it);
            msg = sprintf('Option number: %d \nTime remaining: %4.1fs\n', it, remaining);
            fprintf([reverseStr, msg]);
            reverseStr = repmat(sprintf('\b'), 1, length(msg));
        end
        coords = unique(coords', "rows")';
        geoscatter(coords(1, :), coords(2, :), 'b')
        hold on
%         coords = 90 - coords;
%         geoscatter(coords(1, :), coords(2, :), 'r')
    end