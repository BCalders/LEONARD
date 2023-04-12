function [dopV, fo, r, vel] = dopShift(startTime, stopTime, sat, gs, fe)
% DOPSHIFT - Calculate the Doppler shift between multiple satellites and a stationary ground based observer.
%
% INPUTS: 
%   startTime - StartTime of simulation,
%   stopTime - End time of the simulation,
%   sat - Matrix of satellites to be used in the calculation,
%   gs - Stationary ground based obeserver to be used in the calculation
% 
% OUTPUTS:
%   dopV - Calculates the velocity between SAT and GS,
%   fo - calculates the Doppler shift of the emitted frequency FE,
%   r - returns the range between SAT and GS
% 
% Bas Calders 2022 
% Based on MathWorks source: https://nl.mathworks.com/help/satcom/ug/calculate-latency-and-doppler-in-a-satellite-scenario.html.
% Version 1.1

c = physconst("Lightspeed");                                    % Speed of light

numSats = length(sat);                                          % Number of satellites
simTime = minutes(stopTime - startTime);                        % Simulation time in minutes
ac = access(sat, gs);                                           % Access object
acStatus = accessStatus(ac);                                    % Check if there is a satellite in view of te gs

[az, azall, el, elall] = deal(zeros(simTime+1, numSats));       % Preallocate variables
[satV, satVall, dir] = deal(zeros(3,simTime+1, numSats));
[dopV, r, rall] = deal(nan(simTime+1,numSats));

for iMinute = 0:simTime 
    time = startTime + minutes(iMinute);
    idxMin = iMinute + 1;
    satVis = acStatus(:, idxMin);                               % Check if satellites are visible
    if max(satVis)                                              % only calculating if at least one sat visible
        [azall(idxMin, :),elall(idxMin, :), rall(idxMin, :)] = aer(sat, gs, time);  % Calculate the azimuth, elevation and range
        az(idxMin, logical(satVis)) = azall(idxMin, logical(satVis));               % Only save the values of the visible satellites
        el(idxMin, logical(satVis)) = elall(idxMin, logical(satVis));
        r(idxMin, logical(satVis)) = rall(idxMin, logical(satVis));

        [~,satVall(:, idxMin, :)] = states(sat, time, "CoordinateFrame", "geographic"); % speed in NED(north east down) ~ XYZ
        satV(:, idxMin, logical(satVis)) = satVall(:, idxMin, logical(satVis));     % only save the values of the visible satellites
        dir(:, idxMin, :) = [cosd(el(idxMin, :)).*cosd(az(idxMin, :)); cosd(el(idxMin, :)).*sind(az(idxMin, :)); -sind(el(idxMin, :))]; % direction of the groundstation wrt gs
        dopV(idxMin, :) = dot(satV(:, idxMin, :),dir(:, idxMin, :));                % doppler velocity in m/s
end

dopV(dopV(:, :, :) == 0) = nan;                                 % Set all zeros to nan
fo = (((c ./ (c + dopV)) * fe));                           % Calculate the Doppler shift of the emitted frequency FE
vel = squeeze(satV());                                          % Calculate the velocity between SAT and GS
vel(vel(:, :) == 0) = nan;                                      % Set all zeros to nan
end

