function stimulus=stimcreator(freq,pulselength,totallength)
% Returns n x 2 -matrix that contains time and stimuluscurrent
% freq is frequency (Hz), pulselength is length of pulse (ms),
% totallength is the total length of the stimulus current time course (ms)
% and amplitude is amplitue of the pulse in pA/pF

totallength = totallength*1000; % s to ms conversion
pulselength = pulselength*1000; % s to ms conversion

% WARNING This isn't absolutely accurate stimulator
clear stimulus
endtime=0;
i=2;
amplitude = -2500;

stimulus = zeros(4*ceil(freq)*ceil(totallength)/1000+1,2);

for t=1:totallength
    if mod(t,  round(1000/freq) ) == 0 % Approximation via round has to be done to make sure that mod will work
        starttime = t;
        endtime = t + pulselength;
        stimulus(i,:) = [ starttime-0.0001 0];
        i=i+1;
        stimulus(i,:) = [ starttime amplitude];
        i=i+1;
        stimulus(i,:) = [ endtime amplitude];
        i=i+1;
        stimulus(i,:) = [ endtime+0.0001 0];
        i=i+1;
    end
end

% Correct the time, so that first stimulus is at t=0.00001
% If you put it at t=0 => no AP ???
stimulus(:,1) = stimulus(:,1) - stimulus(3,1);
stimulus = stimulus(3:end,:);
stimulus(1,1) = stimulus(1,1) + 0.00001;

stimulus(:,1) = stimulus(:,1)./1000; %ms to s conversion