

function [X,Y,f,measurements] = Newmark2D_acoustic_fcn(objg,speedmms,defZeroInXsteps,defZeroInYsteps,...
    xmin,xmax,ymin,ymax,dstep,frequencysamples,fstart,fstop,switches,IFbandwidth,calfile,measureSparam)

dbstop if error 
%IF bandwidth is set to
% IFbandwidth

%cal file is set to
% calfile

% %select measurment
% if measureSparam==1
%     PortsToMeas='MeasS31'
% elseif measureSparam==2
%     PortsToMeas='MeasS32'
% elseif measureSparam==3
%     PortsToMeas='MeasS21'
% end

%% initialize Chirp
SamplingRate=44.1E3;
dt=1/SamplingRate;       %44 KHz sampling rate

%chirp
f0    = (7E3);       %start
f1    = (20E3);      %stop
tau   = 0.205;           %chirp ramp time

T_instrDelay=.2;  %delay in start of chirp by instrument
T_startJitter=0.0152; % variation in start or chirp
T_range=0.018; % 3 m range
T_delay=0; %extra delay
T=tau+T_delay+T_range+T_startJitter+T_instrDelay;           %Total Signal time chirp+roundtrip+instrument delay+start_delay
t=(0:dt:T-dt); 
 
signal = zeros(size(t));

%linear chirp function
phi=@(t) 2*pi*(f0*t+(f1-f0)/(2*tau).*t.^2);

chirp = @(t,phi) ((t<T_delay | t>T_delay+tau).*0 ...
                    +((t> T_delay) & (t<=T_delay+tau)).*exp(1j.*phi));

signal=chirp(t,phi(t));

%output signal
%signalOut=[real(chirp_sig.'),zeros(size(t)).']; 
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
fftchirp=1/T*ifftshift(fft(fftshift(signalOut)));
plot(freqs.',real(fftchirp(:,1)));
title('Channel 1: Frequency Domain');
xlim([min(0),max(1/(2*dt))]);

subplot(2,2,4)
freqs=(-1/(2*dt):1/(T):1/(2*dt)-1/(T)).';
fftchirp=1/T*ifftshift(fft(fftshift(signalOut)));
plot(freqs.',real(fftchirp(:,2)));
title('Channel 2: Frequency Domain');
xlim([min(0),max(1/(2*dt))]);

%% initialize instruments


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
s.Rate = SamplingRate;

plotFFT = @(src, event) fftData(event.Data, src.Rate,src.NotifyWhenDataAvailableExceeds, hp_ch1_tdomain, hp_ch1_fdomain, hp_ch2_tdomain, hp_ch2_fdomain);
% test= @(src, event) assignin('base','micData', event.Data);

%collect data and process
% h1 = addlistener(s, 'DataAvailable', test);
h2 = addlistener(s, 'DataAvailable', plotFFT);
% 
figure(hf);



%% setup scan 
nx=length(xmin:dstep:xmax);
ny=length(ymin:dstep:ymax);
clear scanData;
measurements=cell(nx, ny);

[X, Y] = meshgrid(xmin:dstep:xmax,ymin:dstep:ymax);
%[X, Y] = meshgrid(xmin:dstep:xmax,ymax:-dstep:ymin);

stops = size(Y,1)*size(X,2);

%% begin scan!

stopscomp = 0;
remindersent = false;
for yn=1:size(Y,1)
    direction = 2*mod(yn,2)-1;
    if direction>0
        xindex = 1:size(X,2);
    else
        xindex = size(X,2):-1:1;
    end
    for xn=xindex
        tic
        x = X(yn,xn); %no flipdim here--we are interested flipped position anyway
        %y = Y(yn,xn); 
        y = Y(yn,xn); %to account for meshgrid making array that starts at neg pos and goes pos
        Newmark2D_stage_moveToAbsolute(objg,speedmms,defZeroInXsteps,defZeroInYsteps,x,y); %recomm. speed 25 mm/sec

        queueOutputData(s, signalOut);

        s.NotifyWhenDataAvailableExceeds=length(t);
        micData=startForeground(s);
        measurements{yn,xn}=micData;
          
        stopscomp = stopscomp+1;
        timere = (stops-stopscomp)*toc/3600;
        disp(['Est. time remaining: ' num2str(timere) 'hours'])
        
    end
end
%         if timere<=0.75 && ~remindersent
%             remindersent = true;
%             mySMTP = 'smtp.duke.edu';
%             myEmail = {'jjdhunt@gmail.com','jgollub@gmail.com','drsmith@ee.duke.edu'};
%             setpref('Internet','SMTP_Server',mySMTP);
%             setpref('Internet','E_mail',myEmail);
%             recipient = myEmail;
%             subj = 'Near-field scan is nearly complete';
%             msg = ['Hi Guys,',char(10),char(10),'Just a reminder that your current near-field scan will finish in about 45 minutes.',char(10),char(10),'Can''t wait to see you!',char(10),'-Near-field scanner'];
%             sendmail(recipient,subj,msg);
%         end
        tic

%     if switches>6
%         figure(2)
%         imagesc(abs(measurements(:,:,round(length(measurements(1,1,:,1))/2),1)))
%         figure(3)
%         imagesc(abs(measurements(:,:,round(length(measurements(1,1,:,1))/2),7)))
%     elseif (switches<=6) && (switches>0)
%         figure(2)
%         imagesc(abs(measurements(:,:,round(length(measurements(1,1,:,1))/2),1)))
%     else
%         figure(2)
%         imagesc(abs(measurements(:,:,round(length(measurements(1,1,:))/2))))
%     end

%save data
data.X = X;
data.Y = Y;
data.measurements = measurements;

data.chirp.T=T;
data.chirp.dt=dt;
data.chirp.t=t;
data.chirp.f0=f0;
data.chirp.f1=f1;
data.chirp.tau=tau;
data.chirp.signal=signal;
data.chirp.T_instrDelay=T_instrDelay; 
data.chirp.T_startJitter=T_startJitter; 
data.chirp.T_range=T_range;
data.chirp.T_delay=T_delay;

data.samplingRate=SamplingRate;
data.note.order='Data ordered in cells as (yn,xn)';

save(['D:\Acoustic Scanning Data\results-' date],'data','-v7.3')

 delete(h2);
 
 dbclear if error 
