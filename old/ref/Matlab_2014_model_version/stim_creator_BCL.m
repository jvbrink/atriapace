function stimulus = stim_creator_BCL(BCL, pulselength, totallength)
%**************************************************************************
% Returns n x 2 -matrix that contains time and stimuluscurrent.
% BCL is basic cycle length (ms), pulselength is length of pulse (ms),
% totallength is the total length of the stimulus current time course (ms)
% and amplitude is amplitue of the pulse in pA/pF
%
% WARNING This isn't absolutely accurate stimulator
%
%--------------------------------------------------------------------------
% 13.12.2005, Topi Korhonen
% 14.02.2011, Jussi Koivum???ki (changed input from freq to BCL)
%**************************************************************************

totallength = totallength*1000; % s to ms conversion
pulselength = pulselength*1000; % s to ms conversion

clear stimulus
endtime = 0;
i = 2;
amplitude = -1410*2; % twice the threshold

stimulus = zeros(4*ceil(1000/BCL)*ceil(totallength/1000) + 1, 2);

for t = 1:totallength
    if mod(t, round(BCL)) == 0 % Approximation via round has to be done to make sure that mod will work
        starttime = t;
        endtime = t + pulselength;
        stimulus(i,:) = [ starttime-0.0001 0];
        i = i + 1;
        stimulus(i,:) = [ starttime amplitude];
        i = i + 1;
        stimulus(i,:) = [ endtime amplitude];
        i = i + 1;
        stimulus(i,:) = [ endtime+0.0001 0];
        i = i + 1;
    end
end

% Correct the time, so that first stimulus is at t=0.00001
% If you put it at t=0 => no AP ???
stimulus(:,1) = stimulus(:,1) - stimulus(3,1);
stimulus = stimulus(3:end,:);
stimulus(1,1) = stimulus(1,1) + 0.00001;

stimulus(:,1) = stimulus(:,1)./1000; %ms to s conversion