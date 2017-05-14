function fftData(data, Fs, samples,plotHandle_ch1_tdomain, plotHandle_ch1_fdomain,plotHandle_ch2_tdomain, plotHandle_ch2_fdomain)
% Calculate FFT(data) and update plot with it. 

lengthOfData = length(data);
nextPowerOfTwo = 2 ^ nextpow2(lengthOfData); % next closest power of 2 to the length

plotRange = nextPowerOfTwo / 2; % Plot is symmetric about n/2
yDFT = fft(data, nextPowerOfTwo); % Discrete Fourier Transform of data note by not using fftshift we get postive freqs in first half of set

h = yDFT(1:plotRange,:);
abs_h = abs(h);

% assignin('base','dataTEST', data);

freqRange = (0:nextPowerOfTwo-1) * (Fs / nextPowerOfTwo);  % Frequency range
gfreq = freqRange(1:plotRange);  % Only plotting upto n/2 (as other half is the mirror image)
set(plotHandle_ch1_fdomain, 'ydata', db(abs_h(:,1)), 'xdata', gfreq); % Updating the plot
set(plotHandle_ch2_fdomain, 'ydata', db(abs_h(:,2)), 'xdata', gfreq); % Updating the plot

gtime=0:1/Fs:1/Fs*double(samples)-1/Fs;
set(plotHandle_ch1_tdomain, 'ydata', abs(data(:,1)).^2, 'xdata', gtime); % Updating the plot
set(plotHandle_ch2_tdomain, 'ydata', abs(data(:,2)).^2, 'xdata', gtime); % Updating the plot

drawnow; % Update the plot
end


