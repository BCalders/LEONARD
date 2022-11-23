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

c = physconst("Lightspeed");

numSats = length(sat);
simTime = minutes(stopTime - startTime);
ac = access(sat, gs);
acStatus = accessStatus(ac);        % Check if there is a satellite in view of te gs

[az, azall, el, elall] = deal(zeros(simTime+1, numSats));
[satV, satVall, dir] = deal(zeros(3,simTime+1, numSats));
[dopV, r, rall] = deal(nan(simTime+1,numSats));

for iMinute = 0:simTime
    time = startTime + minutes(iMinute);
    idxMin = iMinute + 1;
    satVis = acStatus(:, idxMin);
    if max(satVis)  % only calculating if at least one sat visible
        [azall(idxMin, :),elall(idxMin, :), rall(idxMin, :)] = aer(sat, gs, time);
        az(idxMin, logical(satVis)) = azall(idxMin, logical(satVis));
        el(idxMin, logical(satVis)) = elall(idxMin, logical(satVis));
        r(idxMin, logical(satVis)) = rall(idxMin, logical(satVis));

        [~,satVall(:, idxMin, :)] = states(sat, time, "CoordinateFrame", "geographic"); % ------------speed in NED(north east down) ~ XYZ
        satV(:, idxMin, logical(satVis)) = satVall(:, idxMin, logical(satVis));
        dir(:, idxMin, :) = [cosd(el(idxMin, :)).*cosd(az(idxMin, :)); cosd(el(idxMin, :)).*sind(az(idxMin, :)); -sind(el(idxMin, :))]; % --------- dir of gs irt sat in ned
        dopV(idxMin, :) = dot(satV(:, idxMin, :),dir(:, idxMin, :)); % velocity along the line between gs and sat
    end
end

dopV(dopV(:, :, :) == 0) = nan;
fo = (((c ./ (c + dopV)) * fe) - fe);
vel = squeeze(satV());
vel(vel(:, :) == 0) = nan;
end

