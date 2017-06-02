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
hold on;

sig=chirp(t,phi(t));

plot(t, imag(chirp_r(t-r/c,phi(t), phi_r(t-r/c,r))),'b')
plot(t, imag(chirp_r(t-r/c,phi(t), phi_r(t,r))),'g')
sig_reverse=conj(chirp(t-r/c,phi(t-r/c)));

plot(t, imag(chirp(t-r/c,phi(t))),'b')
plot(t, imag(conj(chirp(t-r/c,phi(t-r/c)))),'r')

plot(t, abs(chirp(t-r/c,phi(t-r/c)).*conj(chirp(t-r/c,phi(t-r/c)))),'k')

figure(102)
       
apollo=fft(real(chirp(t,phi(t))));
plot(0:1/T:44.1E3/2-1/T, apollo(1:end/2),'r')

%% load data and set parameters Updated
c=343; % speed of sound m/s
% 
%  load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
%                   'ScanSpeaker\results-02-May-2017_singleMicOnStage.mat']);
%  load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
%                   'ScanSpeaker\results-15-May-2017SingleMicFixedChirp400x400.mat']);
%  load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
%                   'ScanSpeaker\results-16-May-2017_FixedChirpReally.mat']);
  load(['D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\',...
                   'ScanSpeaker\results-19-May-2017OkActuallyFixedChirp.mat']);             

%chirp Info
T=data.chirp.T;
t=data.chirp.t;
f0=data.chirp.f0;
f1=data.chirp.f1;
% signal=data.chirp.chirp_sig;
signal=data.chirp.signal;
tau=data.chirp.tau;
srate=data.samplingRate;
dt=t(2)-t(1);

%frequencies
freqs=0:1/T:(srate/2)-1/T; % standard srate=44.1KHz
% freqs=0:1/T:(44.1E3/2)-1/T;

phi=@(t) 2*pi*(f0*t+(f1-f0)/(2*tau).*t.^2);

chirp = @(t,phi) ((t<(T-tau)/2 | t>(T+tau)/2).*0 ...
                    +(t> (T-tau)/2 & t<=(T+tau)/2).*exp(1j.*phi));


%plot chirp generated Vs defined chirp
figure(1); clf; 
subplot(1,3,1); hold on; 
plot(t,real(signal));  plot(t,real(chirp(t,phi(t))),'--k'); xlim([0-.1,T+.1]); 
title('Generated Chirp Vs Function')
xlabel('t(s)'); ylabel('Amplitude (a.u.)')
subplot(1,3,2:3); hold on; 
plot(t,real(signal));  plot(t,real(chirp(t,phi(t))),'--k'); 
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
axis tight; xlabel('time (s)'); ylabel('voltage (a.u.)'); xlim([.45,.46]);
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
choosefreq=find(freqs>7000 & freqs<22000);
choosefreq=choosefreq(1:100:end);

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
        measurement_t_offset(ny,nx)=dt*indx; %Changed order (nx, ny)-> (ny, nx) !!!!! but this isn't used again right now
        
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

ChooseFreq=19000;

[~,nn]=min(abs(f-ChooseFreq));
singleFreqData=measurementsfft(:,:,nn);
singleFreq=f(nn);




figure(10); clf;
subplot(2,3,1)
imagesc(X(1,:),Y(:,1),abs(singleFreqData)); axis equal; axis tight;
title('Abs(p)')
subplot(2,3,2)
imagesc(X(1,:),Y(:,1),real(singleFreqData)); axis equal; axis tight;
title('Real(p)')
subplot(2,3,3)
imagesc(X(1,:),Y(:,1),imag(singleFreqData)); axis equal; axis tight;
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
kzs=real(sqrt(k0.^2-kx.^2-ky.^2));

figure(11);
zd=0:.01:1;

% for ii=1:length(zd)
% Exy=(ifft2(fftshift(Ekxky.*exp(1j*kzs*zd(ii))))); %wave is traveling into -z henze positive phase to reverse
% 
% imagesc(abs(Exy))
% axis equal; axis tight;
% title(['xd= ', num2str(zd(ii),'%.3f')])
% drawnow
% pause(.2)
% end

%---------------------------
figure(12); cla; clear mov fig map;
nFrames = length(zd);

Exy=(ifft2(fftshift(Ekxky.*exp(1j*kzs*zd(1))))); 
imagesc(X(1,:),Y(:,1),abs(Exy));

clear mov fig map;
set(gcf, 'Color','white')
set(gca, 'nextplot','replacechildren', 'Visible','off');

fig = getframe(gcf);
[fig,map] = rgb2ind(fig.cdata, 256, 'nodither');
mov = repmat(fig, [1 1 1 nFrames]);

%# create movie
for ii=1:nFrames
Exy=(ifft2(fftshift(Ekxky.*exp(1j*kzs*zd(ii))))); 

imagesc(X(1,:),Y(:,1),abs(Exy));
axis equal; axis tight; 
title(['xd= ', num2str(zd(ii),'%.3f')])
axis on;
drawnow
    fig = getframe(gcf);
    mov(:,:,1,ii) = rgb2ind(fig.cdata, map, 'nodither');
end
close(gcf)

%# create GIF and open
imwrite(mov, map, 'D:\Dropbox (Duke Electric & Comp)\Work\Weekly Slide Data\Acoustics\Movies\BackPropagation.gif', 'DelayTime',0, 'LoopCount',inf)
winopen('D:\Dropbox (Duke Electric & Comp)\Work\Weekly Slide Data\Acoustics\Movies\BackPropagation.gif')


%% Phase reversal technique
%positions
% zrange=0:c/(2*BW):1;

zrange=.1:.1:.6;

[xgrid,ygrid,zgrid]=meshgrid(X(1,:),Y(:,1),zrange);
fest=zeros(size(xgrid));

z0=0;
[x0,y0]=meshgrid(X(1,:),Y(:,1));
% [x0,y0]=meshgrid(X(1,end/4:3*end/4),Y(end/4:3/4*end,1));

clear r;

% tmat=permute(t,[1 3 2]);
% tmat=repmat(tmat,[size(x0),1]);

tIndx=permute(1:length(t),[1 3 2]);

tchirp=chirp(t,phi(t));
tchirp=permute(tchirp,[1 3 2]);
tchirp=repmat(tchirp,size(x0));

shiftMax=round(max(max(sqrt((xgrid(1,1,1)-x0).^2+(ygrid(1,1,1)-y0).^2+(zgrid(1,1,1)-z0).^2)-zgrid(1,1,1)))/c/dt);

delayedChirp=zeros([size(x0),length(t)-shiftMax]);

I=size(xgrid(1,:,1),2);
J=size(xgrid(:,1,1),1);
K=size(xgrid(1,1,:),3);




for i=1:I
    for j=1:1%J
        tic
        for k=1:K
            r=sqrt((xgrid(j,i,k)-x0).^2+(ygrid(j,i,k)-y0).^2+(zgrid(j,i,k)-z0).^2)-zgrid(end/2,end/2,k);
%             r=repmat(r,1,1,length(t));
%         fest(j,i,k)=abs(sum(sum(sum(measurements.*conj( chirp_r(tmat-r/c,phi(tmat-r/c),phi_r(tmat-r/c,r)) )))));
         tIndxShift=tIndx+round((r/c)/dt);
toc
         delayedChirp=tchirp(tIndxShift);
toc     
          fest(j,i,k)=abs(sum(sum(sum(measurements.*conj(delayedChirp)))));
toc          
        end
        timeleft=toc*((I*J)-((i-1)*J+j))/60;
        percent_done=((i-1)*J+j)/(I*J)*100;
    fprintf(['Finished ', num2str(percent_done),'%% \n']);
    fprintf(['Time left ', num2str(timeleft),' mins \n']);
    end
end



    figure(1); subplot(2,1,2); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',fest.^6,'xdata',X(1,:),'Ydata',Y(:,1),'Zdata',zrange);
    axis equal; axis tight; view(3);
    zlabel('downrange (m)'); xlabel('x crossrange (m)');  ylabel('y crossrange (m)');




%% Scan setup
z_p=0; 

%upsample in x and y

pad=2^(nextpow2(max(numel(X(1,:)),numel(Y(:,1)))));
% pad=2^7;

Lx_pad=dx*(pad-1);
Ly_pad=dy*(pad-1);

% calc kx & ky vector information
% [kux,kuy,~]=meshgrid(-(2*pi)/(2*dx):2*pi/(Lx_pad):(2*pi)/(2*dx),...
%                      -(2*pi)/(2*dy):2*pi/(Ly_pad):(2*pi)/(2*dy),...
%                      1:numel(f));

[kux,kuy,~]=meshgrid(linspace(-(2*pi)/(2*dx),(2*pi)/(2*dx),pad), ...
                     linspace(-(2*pi)/(2*dy),(2*pi)/(2*dy),pad), ...
                     linspace(0,f(end),size(measurementsfft,3)));

%calc free space k vector
k=2*pi*f/c;
% k=linspace(0,2*pi*f(end)/c,size(measurementsfft,3));
 k=repmat(permute(k,[1,3,2]),[pad,pad,1]);

%calculate plane wave decomposition (FFT). Note we have to be careful about
%shifting in fft because we are only acting along two of the dimensions.
%Therefor it is easier to use circshift (as opposed to ifftshif which acts
%on all dimentions.

% Sxy=(1/(Lx_pad*Ly_pad))*circshift(fft2(measurements,pad,pad),[pad/2,pad/2]);
Skxky=(1/(Lx_pad*Ly_pad))*circshift(fft2(measurementsfft,pad,pad),[pad/2,pad/2]);

shiftx=dx*(pad/2-ceil(xnum/2));
shifty=dy*(pad/2-ceil(ynum/2));

%padding adds zeros at end of both dimensions, we must phase shift such that panel is centered
Skxky=Skxky.*exp(-1.0j*kux*(shiftx)).*exp(-1.0j*kuy*(shifty));

    %plot decomposition
    figHandle=figure(20);
    scrn = get( groot, 'Screensize'); scrn(1)=2*scrn(3)/3;  scrn(3)=scrn(3)/3;
    set(figHandle,'Position',scrn); clf; subplot(2,2,1);
    imagesc(-Lx_pad/2:dx:Lx_pad/2,-Ly_pad/2:dy:Ly_pad/2,abs(ifft2(ifftshift(Skxky(:,:,60)))));
    title('Padded Input Fields f(1)'); axis equal; axis tight; xlabel('x (m)'); ylabel('y (m)')
    subplot(2,2,2);
    imagesc(kux(1,:),kuy(:,1),abs(Skxky(:,:,60)));
    title('FFT of Fields f(1)'); axis equal; axis tight; xlabel('kx (m)'); ylabel('ky (m)')
    drawnow
    
%calculate min max kz wavenumber
% Kz=2*k;
% kuz=sqrt((2*k).^2-kux.^2-kuy.^2);

%calculate min max kz wavenumber

kuz=real(sqrt((k).^2-(kux).^2-(kuy).^2));
Kz=linspace(min(kuz(:)),max(kuz(:)),size(measurementsfft,3));

% imagesc(X(1,:), linspace(0,f(end),size(measurementsfft,3)),squeeze(abs(Skxky(:,:,300))))

%interpolate to evenly spaced grid
Srmg=zeros(size(kux));
for ii=1:size(kux,2)
    for jj=1:size(kuy,1)
        [~,kRealIndx]=find(kuz(jj,ii,:)~=0);
         if numel(kRealIndx)>2;        
        Srmg(jj,ii,:)=interp1(squeeze(kuz(jj,ii,kRealIndx)),...
            squeeze(Skxky(jj,ii,kRealIndx)),...
            squeeze(Kz),'linear');
         end
    end
end

Srmg(find(isnan(Srmg))) = 0; %set all Nan values to 0

% shift if by location of probe w.r.t. coordinate system
Srmg=Srmg.*exp(1.0j*kuz*z_p);

%apply inverst FFT to get image
fxy = ifftn(fftshift(fftshift(Srmg,1),2));

% fxy = ifftn(fftshift(Srmg.*exp(-1.0j*kuz*z_p),3));
% fxy = ifftn(fftshift(Srmg.*exp(-1.0j*kuz*z_p),3));

image=(abs(fxy)/max(abs(fxy(:))));
 
 %% plot 

% labeling
xx=linspace(-Lx_pad/2,Lx_pad/2,size(fxy,2));
yy=linspace(-Ly_pad/2,Ly_pad/2,size(fxy,1));
zz=linspace(0,c*numel(f)/(BW),numel(f));
% 
% %%. Define spatial domain
% dkz = mean(diff(Kz));
% zz = linspace(0,pi/dkz,size(image,3));
% 
% dkx = mean(diff(linspace(-pi/dx, pi/dx, size(Srmg,1))));
% xx = linspace(-pi/dkx,pi/dkx,size(image,1)); 
% 
% dky = mean(diff(linspace(-pi/dy, pi/dy, size(Srmg,2))));
% yy = linspace(-pi/dky,pi/dky,size(image,2));

%%%
    figure(20); subplot(2,1,2); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',image.^4,'xdata',xx,'Ydata',yy,'Zdata',zz);
    axis equal; axis tight; view(3);
    zlabel('downrange (m)'); xlabel('x crossrange (m)');  ylabel('y crossrange (m)');

 %% plot subimage
xmin=-0.25;
xmax=0.55;
ymin=-0.25;
ymax=0.25;
zmin=0;
zmax=1;

subimage=image(yy>ymin & yy<ymax, xx>xmin & xx<xmax, zz>zmin & zz<zmax);
% subimage=subimage/max(subimage(:));
    figure(20); subplot(2,1,2); cla;
    set(gcf,'color','white'); colormap('parula');
    vol3d('Cdata',subimage.^6,...
        'xdata',xx(xx>xmin & xx<xmax),...
        'Ydata',yy(yy>ymin & yy<ymax),...
        'Zdata',zz(zz>zmin & zz<zmax)...
        );
    zlabel('downrange (m)'); xlabel(' x crossrange (m)');  ylabel('y crossrange (m)');
    axis equal; axis tight; view(3);

