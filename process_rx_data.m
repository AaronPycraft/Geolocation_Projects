cd('C:\Users\Aaron Pycraft\Desktop\EN744_Project_tx_data\adalm_pluto_recordings\');
cd('rx_loc_3')
format compact;
% process rx data
files = dir('*.bb');
fprintf('found %i files\n', size(files));
% fin = files(1).name;
fin = files(1).name;
bbr = comm.BasebandFileReader(fin);
bbr.SamplesPerFrame = 15e6;

data = step(bbr);
release(bbr)
%% get metadata
fc = bbr.CenterFrequency;
fs = bbr.SampleRate;

N = length(data);
t = 1:1:N;
t = t/fs;

f = 1./t;

p = fftshift(fft(data));
spec = 10*log10(abs(p));
% full bandwidth
Low_freq  = (fc-fs/2);
High_freq = (fc+fs/2);
freq = [0:1:N-1]*(fs)/N+Low_freq;
% 10 kHz bw
freq1 = (fc-5e3);
freq2 = (fc+5e3);
zoomIdx = (freq >= freq1) & (freq <= freq2);
freqZoom = freq(zoomIdx);
Nzoom = size(freqZoom);
spectZoom = spec(zoomIdx);


%% Plot time domain
figure(1);
clf(1);

subplot(2,1,1);
plot(t', abs(data))
xlabel('time (sec)');
ylabel('amp');
grid minor;
title('time domain');
ylim([-0.001, max(0.065, 1.1*max(abs(data)))]);

%% Plot freq domain
subplot(2,1,2);
plot(freqZoom/(1e6),spectZoom)
xlabel('freq (MHz)');
ylabel('amp (normalized)')
grid minor;
xlim([freq1, freq2]/1e6)
ylim([min(abs(spectZoom))-10, 10+max(abs(spectZoom))]);
title('spectrum');



%% Find the received time
rise = find(abs(data) > 0.5*max(abs(data)), 1, 'first');             % Index
fall = find((abs(data)) > 0.5*max(abs(data)), 1, 'last');             % Index

data(rise);
time_received = t(rise);
sig_start=t(rise);
sig_end = t(fall);
signal_dur = t(fall)-t(rise);

fprintf('file = %s\n', fin);
fprintf('Signal found @\t%2.3f seconds\n', time_received);
fprintf('Signal start @\t%2.3f seconds\n', sig_start);
fprintf('Signal end @\t%2.3f seconds\n', sig_end);
fprintf('Signal duration\t%2.3f seconds\n', signal_dur);





