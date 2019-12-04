cd('/Users/aaronpycraft/Desktop/adalm_pluto_recordings/');
%% Rx parameters
fc = 437.9999e6
fo = 1e3;
fs = 0.5e6;   % 65105 to 61.4E6 Samp/sec
Gt = -10;
Gr = 60;    

%% Create Radio Objects
sdrdev('Pluto');

deviceSearch = findPlutoRadio
    assert(~isempty(deviceSearch), 'Error, 0 pluto devices found');
if(~exist('rxPluto', 'var'))
rxPluto = sdrrx('Pluto', 'RadioID', deviceSearch(1).RadioID);
    rxPluto.GainSource = 'Manual';
    rxPluto.Gain = Gr;
    rxPluto.CenterFrequency = fc+fo;
    rxPluto.BasebandSampleRate = fs;
    rxPluto.OutputDataType = 'double';  
end

rx_PRI = 15;
tmr = setup_rx_timer(rx_PRI);
tmr.TimerFcn = @(x,y) rx_capture(rxPluto);

% Start recording 2 seconds before expected transmit
nextMin = dateshift(dateshift(datetime('now'), 'start', 'minute', 'next'),'start', 'second', -5)
fprintf('STARTING\n');
fprintf(datestr(nextMin, 'yyyy-mm-DD HH:MM:SS:FFF '));
fprintf('\RX Scheduled\n');

startat(tmr, nextMin);
startTime = tic;    
    
%% Rx Loop

% run once to initialize
% capture(rxPluto, 500);
% rx_capture(rxPluto); 

% Set duration of the transmission PRI repetition
jj=1;
stopRunning = 0; % manually set to end run
try
    while(toc(startTime) < 5*60 && stopRunning==0)
        jj = jj+1;
    end
catch E
    stop(tmr);
    fprintf('\nTimer Stopped\n');    
    % END
    release(rxPluto);
    fprintf(datestr(now, 'yyyy-mm-DD HH:MM:SS:FFF '));
    fprintf('\tEND RECORDING\n');
    rethrow(E);
end

% stop transmission timer
stop(tmr);
fprintf('\nTimer Stopped\n');
        
%% END
release(rxPluto);
fprintf(datestr(now, 'yyyy-mm-DD HH:MM:SS:FFF '));
fprintf('\tEND RECORDING\n');

%% Test Function
%{
% [rxdata,datavalid,overflow] = rxPluto();
rxdata = capture(rxPluto,5,'Seconds','Filename','Recording.bb');

% Close object
pause(1);
release(rxPluto);

%% Process Data
figure(1);
clf(1);

% compute the FFT
N = length(rxdata);
fres = fs/N
spectrum = 10*log(abs(fftshift(fft(rxdata))) / N);   
f = linspace(-fs/2+fres/2, fs/2-fres/2, N); % Create the frequency axis and put the measure in the middle of the bin.
plot(f,spectrum);
xlabel('freq #');
grid minor;

figure(2);
clf(2);

tMad = fs*N
tStart = fs*N
t = 1/fs:1/fs:(1/fs*N);

plot(t, abs(rxdata));
xlabel('sample #');

% bbr = comm.BasebandFileReader('FMRecording.bb');
% bbr.SamplesPerFrame = 4400;
%}

%% Capture & save the spectrum 
function [rxdata] = rx_capture(rxPluto)
    rxTime = now();
    rxTimeStr = datestr(rxTime, 'yyyy-mm-DD HH:MM:SS:FFF ');
    
    fprintf(rxTimeStr);
    fprintf('\tCAPTURING\n');
    fname = 'rxdata';
    captureTime = tic;
    rxdata = capture(rxPluto,2,'Seconds');%,'Filename',fname,'Timestamp', 1);
%     capture(rxPluto,4,'Seconds','Filename',fname, 'Timestamp', 1);
    fprintf(datestr(now(), 'yyyy-mm-DD HH:MM:SS:FFF'));
    fprintf('\tDONE\n');

    pRMS = rms(rxdata)^2;
%     clear rxdata;
    fprintf('\tSAVING\t%s\n', fname, rxTime);
    fprintf('\tPSD = %3.3e\n', pRMS);
    toc(captureTime);
    
%     fprintf('pausing\n');
%     pause();
%     fprintf('resuming\n')
    
    if(pRMS > 0.5e-3)
        
        % Plot data
        fprintf('Plotting\n');
        figure('name', datestr(now()));
        N = length(rxdata);
        fs = rxPluto.BasebandSampleRate;
        fc = rxPluto.CenterFrequency;
        t = [1:1:N]/fs;
        fL = (fc-fs/2);
        fH = (fc+fs/2);
        f = [0:1:N-1] * fs/N + fL;

        subplot(2,1,1);
        plot(t, abs(rxdata));
        title('time domain');
        grid minor; axis tight;

        subplot(2,1,2);
        plot(f, 10*log10(abs(fftshift(fft(rxdata)))))
        title('spectrum');
        grid minor; axis tight;
    else
        fprintf('Power too low. Skipping plot.\n');
    end
    
end

%% Create the rx timer to record the spectrum
function [tmr] = setup_rx_timer(period)
    tmr = timer;
    tmr.Name = 'RX_Timer';
    tmr.BusyMode = 'drop';
    tmr.ExecutionMode = 'fixedDelay';
    tmr.Period = period;
    tmr.StartDelay = 0;
end
