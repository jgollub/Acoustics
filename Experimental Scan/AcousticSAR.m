%% Setup chirp

dt=1/44.1E3;       %44 KHz sampling rate
T=2.2;               %Total Signal time 
t=(0:dt:T-dt);

%chirp
f0    = (7E3);       %start
f1    = (20E3);      %stop
tau   = 2;           %chirp ramp time
t_tau = 0:dt:tau;    
chirp_sig = zeros(size(t));

%linear chirp function
phi=@(t) 2*pi*(f0*t+(f1-f0)/(2*tau).*t.^2);

chirp = @(t,phi) ((t<(T-tau)/2 | t>(T+tau)/2).*0 ...
                    +(t> (T-tau)/2 & t<=(T+tau)/2).*exp(1j.*phi));

signal=chirp(t,phi(t));

signalOut= repmat(real(signal.'), 1, 2);

hf = figure(1);
subplot(2,2,1)
plot(t.',real(signalOut(:,1)));
ylim([-1.1,1.1]);
xlim([min(t),max(t)]);
title('Channel 1: Time domain')

subplot(2,2,2)
plot(t.',real(signalOut(:,2)));
ylim([-1.1,1.1]);
xlim([min(t),max(t)]);
title('Channel 2: Time domain')



subplot(2,2,3)
freqs=(-1/(2*dt):1/(T):1/(2*dt)-1/(T)).';
fftchirp=2/T*ifftshift(fft(fftshift(signalOut)));
plot(freqs.',real(fftchirp(:,1)));
title('Channel 1: Frequency Domain');
xlim([min(0),max(1/(2*dt))]);

subplot(2,2,4)
freqs=(-1/(2*dt):1/(T):1/(2*dt)-1/(T)).';
fftchirp=2/T*ifftshift(fft(fftshift(signalOut)));
plot(freqs.',real(fftchirp(:,2)));
title('Channel 2: Frequency Domain');
xlim([min(0),max(1/(2*dt))]);


%% make and record signal

%channel 1 mic
hf = figure(3);
subplot(2,2,1)
hp_ch1_tdomain = plot(t,zeros(length(t),1));
title(' Channel 1: Time Domain');
xlabel('time (s)')
ylabel('|Y(t)|^2')
grid on;

subplot(2,2,3)
hp_ch1_fdomain = plot(freqs(length(t)/2+1:end),zeros(length(t)/2,1));
title('Channel 1: Discrete FFT Plot');
xlabel('Frequency (Hz)')
ylabel('|Y(f)| dB')
grid on;

subplot(2,2,2)
hp_ch2_tdomain = plot(t,zeros(length(t),1));
title(' Channel 2: Time Domain');
xlabel('time (s)')
ylabel('|Y(t)|^2')
grid on;

subplot(2,2,4)
hp_ch2_fdomain= plot(freqs(length(t)/2+1:end),zeros(length(t)/2,1));
title('Channel 2: Discrete FFT Plot');
xlabel('Frequency (Hz)')
ylabel('|Y(f)| dB')
grid on;

%channel 2 mic

%get devices
d = daq.getDevices;
devMic=d(2); %microphone
devSpeaker=d(4); % speaker

%setup session(sound card as DAQ)
s = daq.createSession('directsound');

% 2 channel output and input sound card
noutchan = 2;
ninchan = 2;

chOut=addAudioOutputChannel(s, devSpeaker.ID, 1:noutchan);
chIn=addAudioInputChannel(s, devMic.ID, 1:ninchan);

%set sampling rate
s.Rate = 44.1E3;

plotFFT = @(src, event) fftData(event.Data, src.Rate,src.NotifyWhenDataAvailableExceeds, hp_ch1_tdomain, hp_ch1_fdomain, hp_ch2_tdomain, hp_ch2_fdomain);
% test= @(src, event) assignin('base','micData', event.Data);

%collect data and process
% h1 = addlistener(s, 'DataAvailable', test);
h2 = addlistener(s, 'DataAvailable', plotFFT);
% 
figure(hf);

nx=2;
ny=1;
clear scanData;
scanData=cell(nx, ny);


for j=1:ny
    for i=1:nx
        
        queueOutputData(s, signalOut);

        s.NotifyWhenDataAvailableExceeds=length(t);
        micData=startForeground(s);
        scanData{i,j}=micData;     
        (j-1)*ny+i
    end
end
%  delete(h1);
 delete(h2);

 figure(4);clf;
 hold on
for i=1:size(scanData,1)
 for j=1:size(scanData,2)
    plot(t, abs(scanData{i,j}(:,1)).^2,'color',rand(1,3))
end
end


%% xcorr calc
figure(5);
v_s=343;
xcorr_val=xcorr(scanData{1}(:,1),scanData{1}(:,2));
t_xcorr=[-flip(t(2:end)),t];
[mp,mv]=max(xcorr_val);
t_shift=t_xcorr(mv);
figure; plot(t_xcorr,xcorr_val);

dist=(t_shift/2)*v_s
(t_shift/2)*v_s/.0254