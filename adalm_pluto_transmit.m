%% TX Loop for gelocation project
% Emits a chirp at a fixed PRI

%% adalm_pluto_transmit
% Transmit test
% Aaron Pycraft
%     clear; clc;
%     format compact;
%% Setup TX
% Tx parameters
fc = 437.9999e6
fs = 2e6;   % 65105 to 61.4E6 Samp/sec
Gt = 0;
Gr = -10;

% Create Radio Objects
sdrdev('Pluto');
deviceSearch = findPlutoRadio;
assert(~isempty(deviceSearch), 'Error, 0 pluto devices found');

txPluto = sdrtx('Pluto', 'RadioID', deviceSearch(1).RadioID);
    txPluto.CenterFrequency = fc;
    txPluto.BasebandSampleRate = fs;
    txPluto.Gain = Gt;

% Setup waveform
tx_waveform = exp(1i*2*pi*fc/fs*(1:50e3)).';
txPluto(tx_waveform)

%% Tx loop

tx_PRI = 15;
tmr = setup_tx_timer(tx_PRI);
tmr.TimerFcn = @(x,y) tx_chirp(fc, fc+5000, fs, txPluto);

nextMin = datestr(dateshift(datetime('now'), 'start', 'minute', 'next'));
fprintf('STARTING\nTX Scheduled to start at %s\n', datestr(nextMin, 'yyyy-mm-DD HH:MM:SS:FFF '));
startat(tmr, nextMin);
startTime = tic;

% Set duration of the transmission PRI repetition
jj=1;
while(toc(startTime) < 60*60)
    jj = jj+1;
end

% stop transmission timer
stop(tmr);
fprintf('\n');
        
%% END
release(txPluto);
fprintf(datestr(now, 'yyyy-mm-DD HH:MM:SS:FFF '));
fprintf('\tEND TRANSMISSION\n');

%% Setup Timer
function [tmr] = setup_tx_timer(period)
    tmr = timer;
    tmr.Name = 'TX_Timer';
%     tmr.stat = 'false';
%     tmr.TimerFcn = @(x,y) tx_chirp
    tmr.BusyMode = 'drop';
    tmr.ExecutionMode = 'fixedDelay';
    tmr.Period = period;
    tmr.StartDelay = 1;
%     start(tmr)
%     startat(now)
end

%% Setup Chirp
% This version creates a 5kHz chirp over 2 seconds
%{
function tx_chirp(freqStart, freqStop, fs, txPluto)
    % Transmit Chirp
    runtime = tic;
    ii=0;
    freq = freqStart;
    
    fprintf(datestr(now, 'yyyy-mm-DD HH:MM:SS:FFF '));
    fprintf('\tTRANSMITTING\n');
    while (toc(runtime) < 2) && (freq < freqStop)
        ii = ii+1;
        % TX
        tx_waveform = exp(1i*2*pi*freq/fs*(1:20e3)).';
        freq = freq + 25;
        txPluto(tx_waveform);
    end
    fprintf('\tfstart=\t%2.2f\n', freqStart);
    fprintf('\tfstop=\t%2.2f\n', freqStop);
%     fprintf('ii=%i\n', ii);

end
%}
function tx_chirp(freqStart, freqStop, fs, txPluto)
    % Transmit Chirp
    runtime = tic;
    ii=0;
    freq = freqStart;
    
    fprintf(datestr(now, 'yyyy-mm-DD HH:MM:SS:FFF '));
    fprintf('\tTRANSMITTING\n');
    while (toc(runtime) < 0.5) && (freq < freqStop)
        ii = ii+1;
        % TX
        tx_waveform = exp(1i*2*pi*freq/fs*(1:50e3)).';
        freq = freq + 25;
        txPluto(tx_waveform);
    end
    fprintf('\tfstart=\t%2.2f\n', freqStart);
    fprintf('\tfstop=\t%2.2f\n', freqStop);
%     fprintf('ii=%i\n', ii);

end
%% Old unused functions
function txTest()
    print('fix me');
    %{
    rxPluto = sdrrx('Pluto', 'RadioID', deviceSearch(2).RadioID,...
        'CenterFrequency', fc, 'BasebandSampleRate', fs,...
        'Gain', Gr ...
        );
    %}

    %% Create Waveform
    %{
    sw = dsp.SineWave;
    sw.Amplitude = 0.5;
    sw.Frequency = 439e6;
    sw.ComplexOutput = true;
    sw.SampleRate = fs;
    sw.SamplesPerFrame  = 5000;
    txWaveform = sw();
    %}
    %{
    chrp = dsp.Chirp;
    chrp.SweepDirection = 'Unidirectional';
    chrp.InitialFrequency = 439e6;
    chrp.TargetFrequency  = 440e6;
    chrp.SweepTime = 1;
    chrp.SampleRate = fs;
    chrp.SamplesPerFrame  = 5000;
    txWaveform = chrp();
    freq = 439e6;
    %}
    % txPluto(txWaveform);
    %% Enable Transmit
    % transmitRepeat(txPluto, txWaveform)
    % pause(5);
end





