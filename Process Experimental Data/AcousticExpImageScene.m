%-------------------------------------------------------------------------------
%SAS - Acoustic Range Migration Algorithm 
%J. Gollub
%-------------------------------------------------------------------------------
%% generate Chirp function


phi=@(t) 2*pi*(f0*t+(f1-f0)/(2*tau).*t.^2);
% 
% chirp = @(t,phi) ((t<(T-tau)/2 | t>(T+tau)/2).*0 ...
%                     +(t> (T-tau)/2 & t<=(T+tau)/2).*exp(1j.*phi));
% 
% phi=@(t,r) 2*pi*(f0*t+(f1-f0)/(2*tau).*t.^2);

% phi_r=@(t,r) -(2*pi/c)*(f0+(f1+f0)/(2*tau).*t).*r;


chirp = @(t,phi) ((t<(T-tau)/2 | t>(T+tau)/2).*0 ...
                    +(t> (T-tau)/2 & t<=(T+tau)/2).*exp(1j.*phi));


% chirp_r = @(t,phi,phi_r) ((t<(T-tau)/2 | t>(T+tau)/2).*0 ...
%                     +(t> (T-tau)/2 & t<=(T+tau)/2).*exp(1j.*phi)).*exp(1j.*phi_r);
%                 
%                 r=6;       
% cla;
% hold on;
% 
% sig=chirp(t,phi(t));
% 
% plot(t, imag(chirp_r(t-r/c,phi(t), phi_r(t-r/c,r))),'b')
% plot(t, imag(chirp_r(t-r/c,phi(t), phi_r(t,r))),'g')
% sig_reverse=conj(chirp(t-r/c,phi(t-r/c)));
% 
% plot(t, imag(chirp(t-r/c,phi(t))),'b')
% plot(t, imag(conj(chirp(t-r/c,phi(t-r/c)))),'r')
% 
% plot(t, abs(chirp(t-r/c,phi(t-r/c)).*conj(chirp(t-r/c,phi(t-r/c)))),'k')
% 
% figure(102)
%        
% apollo=fft(real(chirp(t,phi(t))));
% plot(0:1/T:44.1E3/2-1/T, apollo(1:end/2),'r')

                

%% load data and set parameters
c=343; % speed of sound m/s
T=2.2;
f=0:1/T:44.1E3/2-1/(2*T);
BW=f(end)-f(1);

% %chirp Info
% T=data.chirp.T;
% t=data.chirp.t;
t=0:1/44.1E3:T;

% f0=data.chirp.f0;
% f1=data.chirp.f1;
% t_tau=data.chirp.t_tau;
% chirp_sig=data.chirp.chirp_sig;
% tau=data.chirp.tau;
% srate=data.samplingRate;
% dt=t(2)-t(1);

%frequencies
 freqs=0:1/T:(44.1E3/2)-1/T; % standard srate=44.1KHz

%%
c=343; % speed of sound m/s   

% 
% load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
%                  'SAS of 16oz Can 5-1-2017\results-29-Apr-2017_16OzCan.mat']);              
% 
load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
                 'Duke Target\results-20-May-2017_DUtarget.mat']); 


% load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
%                  'Duke Target\Second Run\results-26-May-2017_CopperCylinder.mat']); 
% 
% load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\Cylinder\',...
%                  'results-04-Jun-2017_Cylinder_ShortChirp.mat']); 
             

                        
%chirp Info
T=data.chirp.T;
t=data.chirp.t;
f0=data.chirp.f0;
f1=data.chirp.f1;
%  signal=data.chirp.chirp_sig;
signal=data.chirp.signal;
tau=data.chirp.tau;
srate=data.samplingRate;
dt=t(2)-t(1);
data.chirp.instrDelay=.2;
% 
% figure(1); subplot(2,1,1); hold on; plot(t,data.measurements{ceil(end/2),ceil(end/2)}(:,1)); drawnow;
%         
% indxKeep=find(t>data.chirp.instrDelay);
% if length(indxKeep)
%     indxKeep=[indxKeep(1)-1,indxKeep];
% end
% 
% data.measurements=cellfun(@(x) x(indxKeep,:),data.measurements,'UniformOutput',false); 
% t=t(indxKeep);
% T=t(end)-t(1);
% 
% plot(t,data.measurements{ceil(end/2),ceil(end/2)}(:,1));

%frequencies
freqs=0:1/T:(srate/2)-1/T; % standard srate=44.1KHz
% freqs=0:1/T:(44.1E3/2)-1/T;

phi=@(t) 2*pi*(f0*t+(f1-f0)/(2*tau).*t.^2);

chirp = @(t,phi) ((t<(T-tau)/2 | t>(T+tau)/2).*0 ...
                    +(t> (T-tau)/2 & t<=(T+tau)/2).*exp(1j.*phi));


%%

choosefreq=find(freqs>10000 & freqs<17591);
choosefreq=choosefreq(1:2:end);

X=-data.X/1000; %!!!!!!!!!!!!note NEGATIVE sign due to moving target instead of target
dx=X(1,2)-X(1,1); %make sure dx is positive
Lx=abs(X(1,end)-X(1,1));
[ynum,xnum]=size(X); 

Y=-data.Y/1000; %!!!!!!!!!!!!note NEGATIVE sign due to moving target instead of target
dy=Y(2,1)-Y(1,1); %make sure d is positive
Ly=abs(Y(end,1)-Y(1,1)); 

measurements=zeros([size(data.measurements),size(data.measurements{1,1},1)]);
measurement_t_offset=zeros(size(data.measurements));
measurementsfft=zeros([size(data.measurements),length(choosefreq)]);

tRef=data.measurements{1,1}(:,2);
temp=0;
figure(60); hold on;
tic
for nx=1:size(data.measurements,2)
    for ny=1:size(data.measurements,1)
        tData=data.measurements{ny,nx};
  
        %determine offset of singal w.r.t set reference measurement
        [XcorrSign,indx]=max(xcorr(tRef,tData(:,2)));
        
%         plot([-flip(t(1:end-1)),t],xcorr(tRef,tData(:,2))); xlim([-.08,.08]);
%         drawnow; pause(.01);
        
        indx=indx-length(t);
        measurements(ny,nx,:)=tData(:,1);
        measurement_t_offset(ny,nx)=dt*indx;
        
        FFTtdata=fft(tData(:,1));
        FFTtdata=FFTtdata(1:end/2);

        for nf=choosefreq
            temp=temp+1;
            %note time shift is phase shift in frequency domain
            measurementsfft(ny,nx,temp)=exp(-1j*2*pi*freqs(nf)*dt*indx)*FFTtdata(nf); %note order needs to match meshgrid
        end
        temp=0;
    end
end
toc


% load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
%                 'SAS of 16oz Can 5-1-2017\results-29-Apr-2017_BkGrnd.mat']);
            
load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
                'Duke Target\results-19-May-2017_background.mat']);

% load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
%                 'Duke Target\Second Run\results-27-May-2017_bkgrnd.mat']);
% 
% load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\Cylinder\',...
%                 'results-04-Jun-2017_BkGrndShortChirp.mat']);

% figure(1); subplot(2,1,2); hold on; plot(data.chirp.t,data.measurements{ceil(end/2),ceil(end/2)}(:,1));
% data.measurements=cellfun(@(x) x(indxKeep,:),data.measurements,'UniformOutput',false); %%%  
% 
% 
% plot(t,data.measurements{ceil(end/2),ceil(end/2)}(:,1));
            
Bkg_measurements=zeros([size(data.measurements),size(data.measurements{1,1},1)]);
Bkg_measurements_t_offset=zeros(size(data.measurements));
backgroundFFT=zeros([size(data.measurements),length(choosefreq)]);

%use tref from before
temp=0;
tic

figure(60)
hold on;
for nx=1:size(data.measurements,2)
    for ny=1:size(data.measurements,1)
        tData=data.measurements{ny,nx};

        %determine offset of singal w.r.t set reference measurement
        [XcorrSign,indx]=max(xcorr(tRef,tData(:,2)));
        
%         plot([-flip(t(1:end-1)),t],xcorr(tRef,tData(:,2))); xlim([-.08,.08]);
%         drawnow; pause(.01);
        indx=indx-length(t);
        Bkg_measurements(ny,nx,:)=tData(:,1);
        Bkg_measurements_t_offset(ny,nx)=dt*indx;
        
        FFTtdata=fft(tData(:,1));
         FFTtdata=FFTtdata(1:end/2);
                 
        for nf=choosefreq
            temp=temp+1;
            %note time shift is phase shift in frequency domain
            backgroundFFT(ny,nx,temp)=exp(-1j*2*pi*freqs(nf)*dt*indx)*FFTtdata(nf); %note order needs to match meshgrid
        end
        temp=0;
    end
end
toc




targetMeasurement=measurementsfft-backgroundFFT;

% apollo=measurementsfft;

f=freqs(choosefreq);
BW=f(end)-f(1);

figure(30);
subplot(2,1,1)

normalizedata=squeeze(mean(abs(targetMeasurement(1,:,:))));
plot(f,normalizedata)
% 
% 
% aveData=mean(normalizedata)
% for nn=1:length(f)
%     if normalizedata(nn)>.2*aveData
% targetMeasurement(:,:,nn)=targetMeasurement(:,:,nn)./normalizedata(nn);
%     end
% end
% 
% subplot(2,1,2)
% plot(f,squeeze(mean(abs(targetMeasurement(1,:,:)))))

% %---------
% ExperimentFFT=cell(size(data.measurements));
% ExperimentFFT=cellfun(@fft,data.measurements, 'UniformOutput', false);
% 
% 
% 
% ExperimentFFT=cellfun(@(x,y) minus(x(:,1),y(:,1)),ExperimentFFT,backgroundFFT,'UniformOutput',false);
% 
% 
% choosefreq=find(freqs>14000 & freqs<18900);
% choosefreq=choosefreq(1:100:end);
% 
% 
% 
% %-----------
% 
% 
% 
% 
% 
% 
% temp=0; 
% for nx=1:size(ExperimentFFT,1)
%     for ny=1:size(ExperimentFFT,2)
%         for nf=33400:100:39600
%             temp=temp+1;
% measurements(ny,nx,temp)=ExperimentFFT{nx,ny}(nf); %not order needs to match meshgrid
%         end
%     temp=0;    
%     end
% end
% f=f(33400:100:39600);
% BW=f(end)-(1);
% 
% %%probe correction???
% 
% figure(1);
% clf; hold on;
% plot(f,db(ExperimentFFT{floor(end/2),floor(end/2)}(1:end/2,1)),'-b');
% plot(f,db(backgroundFFT{floor(end/2),floor(end/2)}(1:end/2,1)), '-r');
% plot(f,db(ExperimentFFT{floor(end/2),floor(end/2)}(1:end/2,1)), '-k');
% 
% plot(f,db(squeeze(squeeze(measurements(floor(end/2),floor(end/2),:)))), '-k');
% 
% ExperimentFFT(:,1)
% % 
% % figure(1);
% % clf; hold on;
% % plot(f,db(sodaCanfft{50,50}(1:end/2,1)),'-b');
% % plot(f,db(backgroundfft{50,50}(1:end/2,1)), '-r');
% % 
% % plot(f,db(sodaCanfft{50,50}(1:end/2,1)...
% %     -backgroundfft{50,50}(1:end/2,1)), '-g');
% % 
% % plot(f,db(signal{50,50}(1:end/2,1)), '-k');
% 
% x_p=0.025;
% y_p=0.1;
% z_p=0.025;

%% load data and set parameters Updated
c=343; % speed of sound m/s
% 
%  load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
%                   'ScanSpeaker\results-02-May-2017_singleMicOnStage.mat']);

 load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
                  'ScanSpeaker\results-15-May-2017SingleMicFixedChirp400x400.mat']);
              

%chirp Info
T=data.chirp.T;
t=data.chirp.t;
f0=data.chirp.f0;
f1=data.chirp.f1;
t_tau=data.chirp.t_tau;
chirp_sig=data.chirp.chirp_sig;
tau=data.chirp.tau;
srate=data.samplingRate;
dt=t(2)-t(1);

%frequencies
freqs=0:1/T:(srate/2)-1/T; % standard srate=44.1KHz
% freqs=0:1/T:(44.1E3/2)-1/T;

%plot chirp generated Vs defined chirp
figure(1); clf; 
subplot(1,3,1); hold on; 
plot(t,real(chirp_sig));  plot(t,real(chirp(t,phi(t))),'--k'); xlim([0-.1,T+.1]); 
title('Generated Chirp Vs Function')
xlabel('t(s)'); ylabel('Amplitude (a.u.)')
subplot(1,3,2:3); hold on; 
plot(t,real(chirp_sig));  plot(t,real(chirp(t,phi(t))),'--k'); 
xlim([T/2-tau/2-0.001,T/2-tau/2+0.001]); xlabel('t(s)'); 
title('Excerpt')

%plot Measured Vs Generated chirp and Soundcard offsets
figure(2); clf;
subplot(3,1,1);
plot(t,db(data.measurements{floor(end/2),floor(end/2)}(:,2)),t,db(data.measurements{floor(end/2),floor(end/2)}(:,1)))
axis tight; xlabel('time (s)'); ylabel('db(voltage)');
title('Signal, time domain (Generated, Measured)');

subplot(3,1,2)
test=fft(data.measurements{floor(end/2),floor(end/2)}(:,1));
plot(freqs,db(test(1:floor(end/2))));
axis tight; xlabel('Freq (Hz)'); ylabel('db(voltage)');
title('Measured Signal in the Frequency Domain');

subplot(3,1,3); cla; hold on;
for i=1:3
    hold on;
    plot(t,squeeze(squeeze(data.measurements{floor(end/2)+i,floor(end/2)+i}(:,2))))
end
axis tight; xlabel('time (s)'); ylabel('voltage (a.u.)'); xlim([.3,.32]);
title('Input signals (Unexpected Offset)');

%load data scan positions
X=data.X/1000;  
dx=X(1,2)-X(1,1);
Lx=X(1,end)-X(1,1);
[ynum,xnum]=size(X);

Y=data.Y/1000;
dy=Y(2,1)-Y(1,1);
Ly=Y(end,1)-Y(1,1);

% Experiment=cell(size(data.measurements));
% Experiment=cellfun(@(x) x(:,1), data.measurements,'UniformOutput', false);
% ExperimentFFT=cell(size(data.measurements));
% signal=cellfun(@(x) minus(x(:,1),y(:,1)),ExperimentFFT,backgroundFFT,'UniformOutput',false);
% ExperimentFFT=cellfun(@fft,data.measurements, 'UniformOutput', false);
% ExperimentFFT=cellfun(@(x) x(1:end/2,1),ExperimentFFT,'UniformOutput',false);

%choose frequencies
choosefreq=find(freqs>17000 & freqs<18900);
choosefreq=choosefreq(1:10:end);

measurements=zeros([size(data.measurements),size(data.measurements{1,1},1)]);
measurement_t_offset=zeros(size(data.measurements));
measurementsfft=zeros([size(data.measurements),length(choosefreq)]);

tRef=data.measurements{1,1}(:,2);
temp=0;
tic
for nx=1:size(data.measurements,1)
    for ny=1:size(data.measurements,2)
        tData=data.measurements{nx,ny};
       
        %determine offset of singal w.r.t set reference measurement
        [XcorrSign,indx]=max(xcorr(tRef,tData(:,2)));
        indx=indx-length(t);
        measurements(ny,nx,:)=tData(:,1);
        measurement_t_offset(ny,nx)=dt*indx; %fixed order
        
        FFTtdata=fft(tData(:,1));
        FFTtdata=FFTtdata(1:end/2);
        for nf=choosefreq
            temp=temp+1;
            %note time shift is phase shift in frequency domain
            measurementsfft(ny,nx,temp)=exp(-1j*2*pi*freqs(nf)*dt*indx)*FFTtdata(nf); %note order needs to match meshgrid
        end
        temp=0;
    end
end
toc

f=freqs(choosefreq);
BW=f(end)-f(1);


%% propagate fields single frequency
nn=ceil(length(f)/1.1);
singleFreqData=targetMeasurement(:,:,nn);
singleFreq=f(nn);

figure(100); clf;
subplot(2,3,1)
imagesc(X(1,:),Y(:,1),abs(singleFreqData)); axis equal; axis tight; axis xy;
title('Abs(p)')
subplot(2,3,2)
imagesc(X(1,:),Y(:,1),real(singleFreqData)); axis equal; axis tight; axis xy;
title('Real(p)')
subplot(2,3,3)
imagesc(X(1,:),Y(:,1),imag(singleFreqData)); axis equal; axis tight;axis xy;
title('Imag(p)')
hold on;

% figure(11)
% for i=1:100
%     if mod(i,2)
%         imagesc(X(1,:),Y(:,1),real(singleFreqData).^2); axis equal; axis tight;
%         pause(.05)
%     elseif ~mod(i,2)
%         imagesc(X(1,:),Y(:,1),imag(singleFreqData).^2); axis equal; axis tight;
%         pause(.05)
%     end
% end


pad=2^nextpow2(max(size(singleFreqData)));
[kx,ky]=meshgrid(linspace(-2*pi/(2*dx),2*pi/(2*dx),pad),linspace(-2*pi/(2*dy),2*pi/(2*dy),pad));

shiftx=dx*(pad/2-ceil(xnum/2));
shifty=dy*(pad/2-ceil(ynum/2));

Ekxky=circshift(fft2(singleFreqData, pad, pad),[pad/2,pad/2]);
Ekxky=Ekxky.*exp(-1.0j*kx*(shiftx)).*exp(-1.0j*ky*(shifty));

subplot(2,3,4)
imagesc(X(1,:),Y(:,1), abs(Ekxky)); axis tight; axis equal; axis xy;
title(['k-space af f= ',num2str(singleFreq)]);
title('K-Space Abs(P)')
subplot(2,3,5)
imagesc(X(1,:),Y(:,1), real(Ekxky)); axis tight; axis equal; axis xy;
title('K-Space Real(P)')

subplot(2,3,6)
imagesc(X(1,:),Y(:,1), imag(Ekxky)); axis tight; axis equal; axis xy;
title('K-Space Imag(P)')

k0=2*pi*singleFreq/c;
kzs=real(sqrt((2*k0).^2-kx.^2-ky.^2));

figure(11);
zd=.1:.005:.8;

for ii=1:length(zd)
Exy=(ifft2(fftshift(Ekxky.*exp(1j*kzs*zd(ii))))); %wave is traveling into -z henze positive phase to reverse

imagesc(X(1,:),Y(:,1),abs(Exy));
axis equal; axis tight; axis xy;
title(['xd= ', num2str(zd(ii),'%.3f')])
drawnow
pause(.05)
end

d=.585;
Exy=(ifft2(fftshift(Ekxky.*exp(1j*kzs*d)))); %wave is traveling into -z henze positive phase to reverse

imagesc(X(1,:),Y(:,1),abs(Exy))
axis equal; axis tight; axis xy;
title(['xd= ', num2str(d,'%.3f')])
drawnow
pause(.2)


%---------------------------
% figure(12); cla; clear mov fig map;
% nFrames = length(zd);
% 
% Exy=(ifft2(fftshift(Ekxky.*exp(1j*kzs*zd(1))))); 
% imagesc(X(1,:),Y(:,1),abs(Exy)); axis xy;
% 
% clear mov fig map;
% set(gcf, 'Color','white')
% set(gca, 'nextplot','replacechildren', 'Visible','off');
% 
% fig = getframe(gcf);
% [fig,map] = rgb2ind(fig.cdata, 256, 'nodither');
% mov = repmat(fig, [1 1 1 nFrames]);
% 
% %# create movie
% for ii=1:nFrames
% Exy=(ifft2(fftshift(Ekxky.*exp(1j*kzs*zd(ii))))); 
% 
% imagesc(X(1,:),Y(:,1),abs(Exy));
% axis equal; axis tight; 
% title(['xd= ', num2str(zd(ii),'%.3f')])
% axis on;
% drawnow
%     fig = getframe(gcf);
%     mov(:,:,1,ii) = rgb2ind(fig.cdata, map, 'nodither');
% end
% close(gcf)

% %# create GIF and open
% imwrite(mov, map, 'D:\Dropbox (Duke Electric & Comp)\Work\Weekly Slide Data\Acoustics\Movies\BackPropagation.gif', 'DelayTime',0, 'LoopCount',inf)
% winopen('D:\Dropbox (Duke Electric & Comp)\Work\Weekly Slide Data\Acoustics\Movies\BackPropagation.gif')

%% Phase reversal technique
%positions
% zrange=0:c/(2*BW):1;

zrange=0:.1:1;

[xgrid,ygrid,zgrid]=meshgrid(X(1,:),Y(:,1),zrange);
fest=zeros(size(xgrid));

z0=0;
[x0,y0]=meshgrid(X(1,:),Y(:,1));

clear r;
tmat=permute(t,[1 3 2]);

I=size(xgrid(1,:,1),2);
J=size(xgrid(:,1,1),1);
K=size(xgrid(1,1,:),3);

for i=1:I
    for j=1:J
        tic
        for k=1:K
            r=sqrt((xgrid(j,i,k)-x0).^2+(ygrid(j,i,k)-y0).^2+(zgrid(j,i,k)-z0).^2);
            r=repmat(r,1,1,length(t));
%         fest(j,i,k)=abs(sum(sum(sum(measurements.*conj( chirp_r(tmat-r/c,phi(tmat-r/c),phi_r(tmat-r/c,r)) )))));
          fest(j,i,k)=abs(sum(sum(sum(measurements.*conj( chirp(tmat-r/c,phi(tmat-r/c)) )))));
        end
        timeleft=toc*((I*J)-((i-1)*J+j))/60;
        percent_done=((i-1)*J+j)/(I*J)*100;
    end
    fprintf(['Finished ', num2str(percent_done),'%% \n']);
    fprintf(['Time left ', num2str(timeleft),' mins \n']);
end



    figure(1); subplot(2,1,2); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',fest.^6,'xdata',X(1,:),'Ydata',Y(:,1),'Zdata',zrange);
    axis equal; axis tight; view(3);
    zlabel('downrange (m)'); xlabel('x crossrange (m)');  ylabel('y crossrange (m)');




%% Range Migration
z_p=0.6; 

% micLoc_x=0;
% micLoc_y=-.03;
% micLoc_z=.04;


% speakerLoc_x=0;
% speakerLoc_y=0.03;
% speakerLoc_z=0.04;

nn=ceil(size(targetMeasurement,3)/2);
%upsample in x and y

pad=2^(nextpow2(max(numel(X(1,:)),numel(Y(:,1)))));
%   pad=2^8;

Lx_pad=abs(dx*(pad-1)); %ensure that it's positive
Ly_pad=abs(dy*(pad-1)); %ensure that it's positive

% calc kx & ky vector information
% [kux,kuy,~]=meshgrid(-(2*pi)/(2*dx):2*pi/(Lx_pad):(2*pi)/(2*dx),...
%                      -(2*pi)/(2*dy):2*pi/(Ly_pad):(2*pi)/(2*dy),...
%                      1:numel(f));

[kux,kuy,~]=meshgrid(linspace(-(2*pi)/(2*dx),(2*pi)/(2*dx),pad), ...
                     linspace(-(2*pi)/(2*dy),(2*pi)/(2*dy),pad), ...
                     linspace(f(1),f(end),size(targetMeasurement,3))); %dummy spacing
 
                 
%calc free space k vector
% k=2*pi*f/c;
% k=linspace(2*pi*f(1)/c,2*pi*f(end)/c,size(measurementsfft,3));
k=linspace(2*pi*f(1)/c,2*pi*f(end)/c,size(targetMeasurement,3));

% k=linspace(0,2*pi/(dt),size(measurementsfft,3));


k=repmat(permute(k,[1,3,2]),[pad,pad,1]);


%calculate plane wave decomposition (FFT). Note we have to be careful about
%shifting in fft because we are only acting along two of the dimensions.
%Therefor it is easier to use circshift (as opposed to ifftshif which acts
%on all dimentions.

% Sxy=(1/(Lx_pad*Ly_pad))*circshift(fft2(measurements,pad,pad),[pad/2,pad/2]);
Skxky=(1/(Lx_pad*Ly_pad))*circshift(fft2(targetMeasurement,pad,pad),[pad/2,pad/2]);


shiftx=dx*(pad/2-ceil(xnum/2));
shifty=dy*(pad/2-ceil(ynum/2));

%padding adds zeros at end of both dimensions, we must phase shift such that panel is centered
Skxky=Skxky.*exp(-1.0j*kux(:,:,:)*(shiftx)).*exp(-1.0j*kuy(:,:,:)*(shifty));

%plot decomposition
magorphase=@(x) abs(x);

figXY=figure(21);
scrn = get( groot, 'Screensize'); scrn(1)=2*scrn(3)/3;  scrn(3)=scrn(3)/3;
set(figXY,'Position',scrn); clf;
% for nn=1:length(f)
    subplot(3,2,1);
    imagesc(-(-Lx_pad/2:abs(dx):Lx_pad/2),-(-Ly_pad/2:abs(dy):Ly_pad/2),magorphase(ifft2(ifftshift(Skxky(:,:,nn)))));
    axis equal; axis tight; xlabel('x(m)'); ylabel('y (m)');axis xy;
    strTitle=sprintf('fields x-y, \n f = %.3g Hz',f(nn));
    title(strTitle);
    drawnow;
    
    subplot(3,2,2);
    imagesc(kux(1,:),kuy(:,1),magorphase(Skxky(:,:,nn)));
    title('K-space'); axis equal; axis tight; xlabel('kx (1/m)'); ylabel('ky (1/m)'), axis xy;
    drawnow; pause(.1);
% end


figZX=figure(20);
scrn = get( groot, 'Screensize'); scrn(1)=scrn(3)/3;  scrn(3)=scrn(3)/3;
set(figZX,'Position',scrn); clf;
% for nn=1:length(f)
    figure(figZX)
    subplot(3,2,1);
    imagesc(f,-(-Ly_pad/2:abs(dy):Ly_pad/2),magorphase(ifft2(squeeze(ifftshift(Skxky(floor(end/2),:,:))))));
    axis tight; xlabel('f?'); ylabel('x(m)');axis xy;
    strTitle=sprintf('fields X-f');
    title(strTitle);
    drawnow;
    
    subplot(3,2,2);
    imagesc(f,kuy(:,1),magorphase(squeeze(Skxky(floor(end/2),:,:))));
    title('K-space'); axis tight; xlabel('w (1/s)'); ylabel('ky (1/m)'), axis xy;
    drawnow; pause(.1);
% end


%calculate min max kz wavenumber
% Kz=2*k;
% kuz=sqrt((2*k).^2-kux.^2-kuy.^2);

%calculate min max kz wavenumber
% Kz=2*k;
kuz=real(sqrt((2*k).^2-(kux).^2-(kuy).^2));
%  Kz=linspace(-max(kuz(:)),max(kuz(:)),size(measurementsfft,3));
 Kz=linspace(min(kuz(:)),max(kuz(:)),size(targetMeasurement,3));

 %after matched filter (phase center correction)
% fshift=(f1-f0)/2;
% Skxky=Skxky.*exp(1.0j*kuz*(z_p));%.*exp(-1.0j*kux*(.2)).*exp(-1.0j*kuy*(.025));
% Skxky=Skxky.*exp(-1.0j*(2*k-kuz)*(z_p));%.*exp(-1.0j*kux*(.2)).*exp(-1.0j*kuy*(.025));
Skxky=Skxky.*exp(-1.0j*(-kuz)*(z_p));%.*exp(-1.0j*kux*(.2)).*exp(-1.0j*kuy*(.025));



figure(21); 
%  for nn=1:length(f)
 subplot(3,2,3)
 imagesc(-(-Lx_pad/2:abs(dx):Lx_pad/2),-(-Ly_pad/2:abs(dy):Ly_pad/2),magorphase(ifft2(ifftshift(Skxky(:,:,nn))))); 
 axis tight; axis equal; axis xy; 
 strTitle=sprintf('fields x-y, post MF \n f = %.3g Hz',f(nn));
 title(strTitle);
 drawnow;
 subplot(3,2,4)
 imagesc(abs(kux(1,:)),abs(kuy(:,1)),magorphase(Skxky(:,:,nn)));  
 title('k-space x-y, post MF'); axis tight; axis equal; axis xy; 
 drawnow; pause(.1);
% end
% for nn=1:length(f)
    figure(figZX)
    subplot(3,2,3);
    imagesc(f,-(-Ly_pad/2:abs(dy):Ly_pad/2),magorphase(ifft2(squeeze(ifftshift(Skxky(floor(end/2),:,:))))));
    axis tight; xlabel('f?'); ylabel('freq');axis xy;
    strTitle=sprintf('fields X-f after MF');
    title(strTitle);
    drawnow;
    
    subplot(3,2,4);
    imagesc(f,kuy(:,1),magorphase(squeeze(Skxky(floor(end/2),:,:))));
    title('K-space after MF'); axis tight; xlabel('w (1/s)'); ylabel('ky (1/m)'), axis xy;
    drawnow; pause(.1);
% end
 
 %   figure; imagesc(real(kuz(:,:,70)));

%interpolate to evenly spaced grid
% Srmg=zeros(size(kux));
% for ii=1:size(kux,2)
%     for jj=1:size(kuy,1)
%         Srmg(jj,ii,:)=interp1(squeeze(kuz(jj,ii,:)),...
%             squeeze(Sxy(jj,ii,:)),...
%             squeeze(Kz(jj,ii,:)),'linear');
%     end
% end
% 
Srmg=zeros(size(kux));
clear keep
for ii=1:size(kux,2)
    for jj=1:size(kuy,1)
        kRealIndx=find(kuz(jj,ii,:));
        if numel(kRealIndx)>2
                        
            StoltInterp=interp1(squeeze(kuz(jj,ii,kRealIndx)),... %calculate Kz values on a regular grid
                squeeze(Skxky(jj,ii,kRealIndx)),... 
                squeeze(Kz),'nearest');
            
             keep=~isnan(StoltInterp);
% % 
%              Srmg(jj,ii,keep)= ...
%                 squeeze(2*pi./(1j*Kz(keep)))... %matched Filter Term
%                 .*StoltInterp(keep);
          
                       Srmg(jj,ii,keep)=StoltInterp(keep); %W/O MF Term
                       
%                        figure(22); 
%                       plot(squeeze(kuz(jj,ii,kRealIndx)),abs(squeeze(Skxky(jj,ii,kRealIndx))), 'k',...
%                           Kz,squeeze(abs(Srmg(jj,ii,:))),'--r') 
        end
    end
end


% Srmg(find(isnan(Srmg))) = 0; %set all Nan values to 0

% figure(21);
% %    for nn=1:length(Kz)
%     subplot(3,2,5)
%     zpos=linspace(0,2*pi/(Kz(2)-Kz(1)),length(Kz));
%     imagesc(-(-Lx_pad/2:abs(dx):Lx_pad/2),-(-Ly_pad/2:abs(dy):Ly_pad/2),magorphase(ifft2(ifftshift(Srmg(:,:,nn)))));
%     axis tight; axis equal; axis xy;
%     strTitle=sprintf('fields x-y, post MF, Stolt Interp \n z = %.3g (m) ',zpos(nn));
%     title(strTitle);
%     drawnow;
%     subplot(3,2,6)
%     imagesc(abs(kux(1,:)),abs(kuy(:,1)),magorphase(Srmg(:,:,nn)));
%     title('k-space x-y, post MF'); axis tight; axis equal; axis xy;
%     drawnow; pause(.1);
% %   end
% % for nn=1:length(f)
%     figure(figZX)
%     subplot(3,2,5);
%     imagesc(linspace(0,c/(2*(f(2)-f(1))),length(Kz)),-(-Ly_pad/2:abs(dy):Ly_pad/2),magorphase(ifft2(squeeze(ifftshift(Srmg(floor(end/2),:,:))))));
%     axis tight; xlabel('x (m)'); ylabel('freq');axis xy;
%     strTitle=sprintf('fields X-f after MF, after Stolt Interp');
%     title(strTitle);
%     drawnow;
%     
%     subplot(3,2,6);
%     imagesc(Kz,kuy(:,1),magorphase(squeeze(Srmg(floor(end/2),:,:))));
%     title('K-space after MF, after Stolt Interp'); axis tight; xlabel('kz (1/m)'); ylabel('ky (1/m)'), axis xy;
%     drawnow; pause(.1);
% % end

%apply inverst FFT to get image
%   fxy = (ifftn(fftshift(fftshift(fftshift(Srmg,1),2),3)));
%     fxy = ifftshift(ifftn(fftshift(Srmg)));
%         fxy = fftshift(ifftn((Srmg)));

 fxy = (2*BW/c)*(Lx_pad*Ly_pad)*ifftshift(ifftn(Srmg,2*size(Srmg)),3);

% fxy = ifftn(fftshift(fftshift(fftshift(fftshift(Srmg,3),1),2),3));
% fxy = (2*BW/c)*(Lx_pad*Ly_pad)*ifftshift(ifftn(Srmg),3);

% fxy = (2*BW/c)*(Lx_pad*Ly_pad)*fftshift(ifftn(Srmg),3);

% fxy = fftshift(ifftn(fftshift(fftshift(Srmg,1),2)),3);
image=(abs(fxy)/max(abs(fxy(:))));
 
% test=(abs(fxy)<.85).*abs(fxy);
% 
% image=test/max(test(:));

 %% plot 

% %labeling
% xx=linspace(-Lx_pad/2,Lx_pad/2,size(fxy,2));
% yy=linspace(-Ly_pad/2,Ly_pad/2,size(fxy,1));
% zz=linspace(-c*length(Kz)/(2*2*BW),c*length(Kz)/(2*2*BW),length(Kz))+z_p;

% xx=linspace(-Lx_pad/2,Lx_pad/2,size(fxy,2));
% yy=linspace(-Ly_pad/2,Ly_pad/2,size(fxy,1));
% zz=linspace(0,c*length(Kz)/(2*BW),size(fxy,3))+z_p;


dkz = mean(diff(Kz));
xx=linspace(-Lx_pad/2,Lx_pad/2,size(fxy,2));
yy=linspace(-Ly_pad/2,Ly_pad/2,size(fxy,1));
zz=linspace(-pi/dkz,pi/dkz,size(fxy,3))+z_p;

%%. Define spatial domain
% dkz = mean(diff(Kz));
% zz = linspace(0,pi/dkz,size(image,3));
% 
% dkx = mean(diff(linspace(-pi/dx, pi/dx, size(Srmg,1))));
% xx = linspace(-pi/dkx,pi/dkx,size(image,1));
% 
% dky = mean(diff(linspace(-pi/dy, pi/dy, size(Srmg,2))));
% yy = linspace(-pi/dky,pi/dky,size(image,2));


% dkz = mean(diff(Kz));
% zz = linspace(0,pi/dkz,size(image,3));
% 
% dkx = mean(diff(linspace(-pi/abs(dx), pi/abs(dx), size(Srmg,1))));
% xx = linspace(pi/dkx,-pi/dkx,size(image,1));
% 
% dky = mean(diff(linspace(-pi/abs(dy), pi/abs(dy), size(Srmg,2))));
% yy = linspace(pi/dky,-pi/dky,size(image,2));

%%%
    figure(41);% subplot(2,1,2); 
    cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',image.^3,'xdata',xx,'Ydata',yy,'Zdata',zz);
    axis equal; axis tight; view(3);
    zlabel('downrange (m)'); xlabel('x crossrange (m)');  ylabel('y crossrange (m)');
    
    
%     
%         figure(20); subplot(2,1,2); cla;
%     set(gcf,'color','white'); colormap('parula');
%     vol3d('Cdata',image.^5);
%     axis equal; axis tight; view(3);
%     zlabel('downrange (m)'); xlabel('x crossrange (m)');  ylabel('y crossrange (m)');

 %% plot subimage
xmin=min(xx);
xmax=max(xx);
ymin=min(yy);
ymax=max(yy);
 zmin=0.5;
 zmax=.8;

subimage=image(yy>=ymin & yy<=ymax, xx>=xmin & xx<=xmax, zz>=zmin & zz<=zmax);
%  subimage=subimage/max(subimage(:));
    figure(20); subplot(2,1,2); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',subimage.^6.5,...
        'xdata',xx(xx>=xmin & xx<=xmax),...
        'Ydata',yy(yy>=ymin & yy<=ymax),...
        'Zdata',zz(zz>=zmin & zz<=zmax)...
        );
    zlabel('downrange (m)'); xlabel(' x crossrange (m)');  ylabel('y crossrange (m)');
    axis equal; axis tight; view(3);

