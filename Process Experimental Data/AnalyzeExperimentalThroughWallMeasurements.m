
%% load data
path='D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\RFAcousticModulation\NFS-27-Oct-2017_wall_Pole_in_front.mat';
WithAcousticMod=loadAcousticModData(path);

path='D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\RFAcousticModulation\NFS-07-Dec-2017_NoAcousticMod_FixedTrigger.mat';
WithoutAcousticMod = loadAcousticModData(path);

c=3e8;


%% Raw Fields at aperture

%modulated

nn=1
fig_ExpTotalRFwa=figure(20);
plotSelectFields(fig_ExpTotalRFwa, 'up','fields',WithAcousticMod,'measurements',nn);
WithAcousticMod = kComponents(WithAcousticMod,'measurements',nn);
plotSelectFields(fig_ExpTotalRFwa, 'up','kspace',WithAcousticMod,'measurements',nn);
%unmodulated
fig_ExpTotalRFwoa=figure(21);
plotSelectFields(fig_ExpTotalRFwoa, 'down','fields',WithoutAcousticMod,'measurements',nn);
WithoutAcousticMod = kComponents(WithoutAcousticMod,'measurements',nn);
plotSelectFields(fig_ExpTotalRFwoa, 'down','kspace',WithoutAcousticMod,'measurements',nn);

%% propagated fields

%modulated
nn=1;
fig_propagated=figure(22);
WithAcousticMod.d=0.21;
WithAcousticMod = kComponents(WithAcousticMod,'measurements',nn);
WithAcousticMod.k0=2*pi*WithAcousticMod.f(nn)/c;

WithAcousticMod.kzs=real(sqrt((2*WithAcousticMod.k0).^2-WithAcousticMod.kx.^2-WithAcousticMod.ky.^2));
WithAcousticMod.PropagatedExy=(ifft2(fftshift(WithAcousticMod.Ekxky_measurements.*exp(1j*WithAcousticMod.kzs*WithAcousticMod.d)))); %wave is traveling into -z henze positive phase to reverse

plotPropgatedFields(fig_propagated,'up',WithAcousticMod,'PropagatedExy','plottype',@abs)


%unmodulated
nn=1;
fig_propagated=figure(22);
WithoutAcousticMod.d=0.21;
WithoutAcousticMod = kComponents(WithoutAcousticMod,'measurements',nn);
WithoutAcousticMod.k0=2*pi*WithoutAcousticMod.f(nn)/c;

WithoutAcousticMod.kzs=real(sqrt((2*WithoutAcousticMod.k0).^2-WithoutAcousticMod.kx.^2-WithoutAcousticMod.ky.^2));
WithoutAcousticMod.PropagatedExy=(ifft2(fftshift(WithoutAcousticMod.Ekxky_measurements.*exp(1j*WithoutAcousticMod.kzs*WithoutAcousticMod.d)))); %wave is traveling into -z henze positive phase to reverse

plotPropgatedFields(fig_propagated,'down',WithoutAcousticMod,'PropagatedExy','plottype',@abs)


%%  compare modulated and unmodulated

% fieldDifference=WithAcousticMod;
% nn=1;
% d=0.35;
% for nn=1:length(WithAcousticMod.f)
% fieldDifference.measurements=1*WithAcousticMod.measurements(:,:,nn)-1*WithoutAcousticMod.measurements(:,:,nn);
% fieldDifference=kComponents(fieldDifference);
% 
% k0=2*pi*fieldDifference.f/c;
% kzs=real(sqrt((2*k0).^2-fieldDifference.kx.^2-fieldDifference.ky.^2));
% fieldDifference.PropagatedExy = ifft2(fftshift(fieldDifference.Ekxky.*exp(1j*kzs*d)));
% fig_plotDiff = figure(10);
% 
% imagesc(X(1,:),Y(:,1),db(fieldDifference.PropagatedExy),[-60,-40]);
% axis equal; axis tight; axis xy;
% title(['xd= ', num2str(d,'%.3f')])
% set(gca,'XDir','Reverse')
% drawnow; disp(['step', num2str(nn)])
% 
% pause(.1)
% end


%% movie??
% 
% fig_propagated=figure(21)
% WithAcousticMod.d=0.33;
% WithAcousticMod.PropagatedExy=(ifft2(fftshift(WithAcousticMod.Ekxky_measurements.*exp(1j*kzs*WithAcousticMod.d)))); %wave is traveling into -z henze positive phase to reverse
% 
% plotPropgatedFields(fig_propagated,'up',WithAcousticMod,'PropagatedExy','plotType',@db)

% figure()
% 
% WithAcousticMod=kComponents(WithAcousticMod);
% 
% k0=2*pi*WithAcousticMod.f(nn)/c;
% kzs=real(sqrt((2*k0).^2-WithAcousticMod.kx.^2-WithAcousticMod.ky.^2));
% 
% plotSelectFields(fig_ExtractedRF, 'down' ,'kspace', WithAcousticMod,'measurements',1);
% 
% fig_propagated=figure();
% WithAcousticMod.d=.35;
% WithAcousticMod.PropagatedExy=(ifft2(fftshift(WithAcousticMod.Ekxky.*exp(1j*kzs*WithAcousticMod.d)))); %wave is traveling into -z henze positive phase to reverse
% 
% plotPropgatedFields(fig_propagated,'down',WithAcousticMod)



%% characteristic time/freq domain

%modulated
figure(100);

n_indx=9;
subplot(2,1,1);
plot(WithAcousticMod.tMeas,real(squeeze(squeeze(WithAcousticMod.measurements(n_indx,n_indx,:)))));
hold on;
subplot(2,1,2);
fft_one_position=fft(squeeze(WithAcousticMod.measurements(n_indx,n_indx,:)));
plot(WithAcousticMod.CWSamplingFreqs(1:ceil(end/2)),db(fft_one_position(1:ceil(end/2))));
hold on;

%unmodulated
figure(102);
subplot(2,1,1);
plot(WithoutAcousticMod.tMeas,real(squeeze(squeeze(WithoutAcousticMod.measurements(n_indx,n_indx,:)))));
hold on;
subplot(2,1,2);
fft_one_position=fft(squeeze(WithoutAcousticMod.measurements(n_indx,n_indx,:)));
plot(WithoutAcousticMod.CWSamplingFreqs(1:ceil(end/2)),db(fft_one_position(1:ceil(end/2))));
hold on;

%% 
%  delta_f_a=5;
% % gaussFreqFunctional=(1/(delta_f_a*sqrt(2*pi)))*exp(-(freqs-f_a).^2)/(2*delta_f_a.^2);
% 
% BW=20; %Hz
% f_a=WithAcousticMod.acousticSignal.f0;
% 
% sincFilter=2*BW*sinc(2*BW*WithoutAcousticMod.tMeas);
% SqFilter=rectangularPulse(f_a-BW/2,f_a+BW/2,WithoutAcousticMod.measSampleFreqs);
% 
% figure(110)
% subplot(2,1,1)
% plot(WithoutAcousticMod.measSampleFreqs,db(SqFilter));
% 
% subplot(2,1,2)
% sincFilter=2*BW*sinc(2*BW*WithoutAcousticMod.tMeas);
% fftSqWave=ifft(SqFilter);
% 
% plot(WithoutAcousticMod.tMeas,real(fftSqWave),WithoutAcousticMod.tMeas,sincFilter)
%% Extract modulated time domain fields
selectFreq=0;
[~,nn] = min(abs(WithAcousticMod.CWSamplingFreqs-selectFreq))

WithAcousticMod.measurementsfftTimeDomain = fft(WithAcousticMod.measurements,[],3);
figure();
plot(WithAcousticMod.CWSamplingFreqs(1:ceil(end/2)),squeeze(squeeze(db(WithAcousticMod.measurementsfftTimeDomain(40,40,1:ceil(end/2))))))
% plot(squeeze(squeeze(db(WithAcousticMod.measurementsfftTimeDomain(40,40,:)))))
 test=figure();

 
 
WithAcousticMod.d=.21;
WithAcousticMod = kComponents(WithAcousticMod,'measurements',nn);
WithAcousticMod.k0=2*pi*WithAcousticMod.f(nn)/c;
WithAcousticMod.kzs=real(sqrt((2*WithAcousticMod.k0).^2-WithAcousticMod.kx.^2-WithAcousticMod.ky.^2)); 

 frame=struct('cdata',[],'colormap',[]);
 dummy=0;
 for nn=1:250
     dummy=dummy+1;
     imagesc(real(squeeze(WithAcousticMod.measurementsfftTimeDomain(:,:,nn))))
     WithAcousticMod = kComponents(WithAcousticMod,'measurementsfftTimeDomain',nn);
     
%      abs(WithAcousticMod.FreqSlicePropagatedExy)
     WithAcousticMod.FreqSlicePropagatedExy=(ifft2(fftshift(WithAcousticMod.Ekxky_measurementsfftTimeDomain.*exp(1j*WithAcousticMod.kzs* WithAcousticMod.d))));
     plotPropgatedFields(test,'down',WithAcousticMod,'FreqSlicePropagatedExy','clim',[-50,60]);
     title(['Frequency = ',num2str(WithAcousticMod.CWSamplingFreqs(nn)), ' Hz']);
     drawnow;
     
     frame(dummy)=getframe(gcf);
     
     if nn==1
         for tt=1:60
             dummy=dummy+1;
             frame(dummy)=getframe(gcf);
         end
     elseif nn==159

         for tt=1:60
             dummy=dummy+1;
             frame(dummy)=getframe(gcf);
         end
     end
 end
 
 
  video = VideoWriter('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\Presentations\TemporalFrequency.avi'); %video object
 open(video);
 writeVideo(video,frame);
 close(video);
 
 
 %% plot side-by-side

 WithAcousticMod.d=.31;
 fig_sidebyside=figure();
 title(['Frequency = ',num2str(WithAcousticMod.CWSamplingFreqs(nn)), ' Hz']);
      
%  nn=1:size(WithAcousticMod.measurementsfftTimeDomain,3); 
  nn=1;
     WithAcousticMod = kComponents(WithAcousticMod,'measurementsfftTimeDomain',nn);
     WithAcousticMod.FreqSlicePropagatedExy=(ifft2(fftshift(WithAcousticMod.Ekxky_measurementsfftTimeDomain.*exp(1j*WithAcousticMod.kzs* WithAcousticMod.d))));
%      plotPropgatedFields(fig_sidebyside,'sidebyside_left',WithAcousticMod,'FreqSlicePropagatedExy','clim',[-50,60],'plotType',@angle);
      plotPropgatedFields(fig_sidebyside,'sidebyside_left',WithAcousticMod,'FreqSlicePropagatedExy','plotType',@real);

     nn=159:159; %159
     WithAcousticMod = kComponents(WithAcousticMod,'measurementsfftTimeDomain',nn);
     WithAcousticMod.FreqSlicePropagatedExy=(ifft2(fftshift(WithAcousticMod.Ekxky_measurementsfftTimeDomain.*exp(1j*WithAcousticMod.kzs* WithAcousticMod.d))));
%      plotPropgatedFields(fig_sidebyside,'sidebyside_right',WithAcousticMod,'FreqSlicePropagatedExy','clim',[-50,60],'plotType',@angle);
     plotPropgatedFields(fig_sidebyside,'sidebyside_right',WithAcousticMod,'FreqSlicePropagatedExy','plotType',@real);
     
     drawnow;
%%
figurewithout=figure(110);
  
      WithoutAcousticMod.measurementsfftTimeDomain = fft(WithoutAcousticMod.measurements,[],3);
          
     WithoutAcousticMod.d=.21;
     
     WithoutAcousticMod = kComponents(WithoutAcousticMod,'measurementsfftTimeDomain',1);
     WithoutAcousticMod.FreqSlicePropagatedExyDC=(ifft2(fftshift(WithoutAcousticMod.Ekxky_measurementsfftTimeDomain.*exp(1j*WithoutAcousticMod.kzs* WithoutAcousticMod.d))));
%      WithoutAcousticMod.FreqSlicePropagatedExyDC(40:80,50:90)=0;
    plotPropgatedFields(figurewithout,'sidebyside_left',WithoutAcousticMod,'FreqSlicePropagatedExyDC','plottype',@real);
     
     WithoutAcousticMod = kComponents(WithoutAcousticMod,'measurementsfftTimeDomain', 159);
     WithoutAcousticMod.FreqSlicePropagatedExy200Hz=(ifft2(fftshift(WithoutAcousticMod.Ekxky_measurementsfftTimeDomain.*exp(1j*WithoutAcousticMod.kzs* WithoutAcousticMod.d))));
%      WithoutAcousticMod.FreqSlicePropagatedExy200Hz(40:80,50:90)=0;
     plotPropgatedFields(figurewithout,'sidebyside_right',WithoutAcousticMod,'FreqSlicePropagatedExy200Hz','plottype',@real);
     
     
     %% compare DC and AC components 

     
     
     WithAcousticMod.k0=2*pi*WithAcousticMod.f(nn)/c;
WithAcousticMod.kzs=real(sqrt((2*WithAcousticMod.k0).^2-WithAcousticMod.kx.^2-WithAcousticMod.ky.^2)); 

figureDCvs200Hz=figure(102);
nn=1 %DC
     WithAcousticMod = kComponents(WithAcousticMod,'measurementsfftTimeDomain',1);
     WithAcousticMod.FreqSlicePropagatedExyDC=(ifft2(fftshift(WithAcousticMod.Ekxky_measurementsfftTimeDomain.*exp(1j*WithAcousticMod.kzs* WithAcousticMod.d))));
%      WithAcousticMod.FreqSlicePropagatedExyDC(40:80,50:90)=0;
    plotPropgatedFields(figureDCvs200Hz,'sidebyside_left',WithAcousticMod,'FreqSlicePropagatedExyDC','plottype',@db);

    selectFreq=200;
[~,nn] = min(abs(WithAcousticMod.CWSamplingFreqs-selectFreq)); 
     WithAcousticMod = kComponents(WithAcousticMod,'measurementsfftTimeDomain', nn);
     WithAcousticMod.FreqSlicePropagatedExy200Hz=(ifft2(fftshift(WithAcousticMod.Ekxky_measurementsfftTimeDomain.*exp(1j*WithAcousticMod.kzs* WithAcousticMod.d))));
%      WithAcousticMod.FreqSlicePropagatedExy200Hz(40:80,50:90)=0;
     plotPropgatedFields(figureDCvs200Hz,'sidebyside_right',WithAcousticMod,'FreqSlicePropagatedExy200Hz','plottype',@db);
  
      figuresubtract=figure(104);
%      subtract = @(a,b) a-b;
     subproj = @(a,b) a-b*sum(a(:).*conj(b(:)))/sum(b(:).*conj(b(:)));
   
%      WithAcousticMod.DCMinusAC = subproj(WithAcousticMod.FreqSlicePropagatedExyDC,WithAcousticMod.FreqSlicePropagatedExy200Hz);
     
     
%       plotPropgatedFields(figuresubtract,'sidebyside_left',WithAcousticMod,'DCMinusAC','plottype',@real);
    
      
        figuresubtract=figure(105);
      
      plotPropgatedFields(figuresubtract,'sidebyside_left',WithAcousticMod,'FreqSlicePropagatedExyDC','plottype',@db,'clim',[-80,50]);
%      for mult=1:400
%calculate ratio of DC/ACmodulated signal ~magnitude of signal

WithAcousticMod.SubtractedProjection = subproj(WithAcousticMod.FreqSlicePropagatedExyDC,WithAcousticMod.FreqSlicePropagatedExy200Hz);

% WithAcousticMod.threshold = abs(WithAcousticMod.FreqSlicePropagatedExyDC)./abs(WithAcousticMod.FreqSlicePropagatedExy200Hz);
  WithAcousticMod.thresholdComplex =  WithAcousticMod.FreqSlicePropagatedExyDC./WithAcousticMod.FreqSlicePropagatedExy200Hz;
%   WithAcousticMod.ABStar =  (WithAcousticMod.FreqSlicePropagatedExyDC.*conj(WithAcousticMod.FreqSlicePropagatedExy200Hz));
%                          
%   
  
  
     threshold=100;
%      %set all values greater than 100 to zero
% 
%      WithAcousticMod.threshold(WithAcousticMod.weights>threshold)=0;
%        WithAcousticMod.ABStar(WithAcousticMod.weights>threshold)=0;
      
%       plotPropgatedFields(figuresubtract,'sidebyside_right',WithAcousticMod,'ABStar','plottype',@angle,'clim', [-100,100]);
      plotPropgatedFields(figuresubtract,'sidebyside_right',WithAcousticMod,'SubtractedProjection','plottype',@db,'clim',[-80,50]);

      
      drawnow
%      pause
%      end

%% Remove signal from wall using acoustic modulation

% dt = t(2)-t(1);
% T  = WithAcousticMod.sweepTime;


if isfield(WithAcousticMod,'acousticSignal')
    f_a          = WithAcousticMod.acousticSignal.f0;
    delta_f_a    = 20; %Hz
else
    f_a          = 0.2E3; %!!!!!!!!!!!!!!!!! HARD CODED 
end
% f_modulation = @(t) (1i*2*pi*(f_a +t.^2.*(delta_f_a^2/(1i*2))));
BW=5;
B_H = (2*f_a+BW/2);
B_L = (2*f_a-BW/2);
f_modulation = @(t) (2*B_H*sinc(2*B_H*t)- 2*B_L*sinc(2*B_L*t));
figure(102)
plot(WithAcousticMod.CWSamplingFreqs(1:ceil(end/2)),real(fft(f_modulation(WithAcousticMod.tMeas(1:ceil(end/2))))));

Aprime=struct();
Aprime.X = WithAcousticMod.X;
Aprime.Y = WithAcousticMod.Y;
Aprime.f=WithAcousticMod.f;

% Aprime.measurements=removeModulatedRFSignal(WithAcousticMod.measurements,f_modulation,t,WithAcousticMod.X,WithAcousticMod.Y);
  Aprime.measurements=removeModulatedRFSignal(WithAcousticMod.measurements,f_modulation,WithAcousticMod.tMeas,WithAcousticMod.X,WithAcousticMod.Y);

%% Extract Target fields
fig_ExtractedRF=figure(22);
plotSelectFields(fig_ExtractedRF, 'down', 'fields',Aprime,'measurements',1); %!!!!!!!!!!! need to fix

Aprime=kComponents(Aprime,'measurements',159);

k0=2*pi*Aprime.f(nn)/c;
kzs=real(sqrt((2*k0).^2-Aprime.kx.^2-Aprime.ky.^2));

plotSelectFields(fig_ExtractedRF, 'down' ,'kspace', Aprime,'measurements',1);

fig_propagated=figure();

Aprime.d=.31;
Aprime.Exy_Propagated=(ifft2(fftshift(Aprime.Ekxky_measurements.*exp(1j*kzs*Aprime.d)))); %wave is traveling into -z henze positive phase to reverse

plotPropgatedFields(fig_propagated,'down',Aprime,'Exy_Propagated','plotType',@db);

% imagesc(X(1,:),Y(:,1),abs(Aprime.PropagatedAxy))
% axis equal; axis tight; axis xy;
% title(['xd= ', num2str(d,'%.3f')])
% set(gca,'XDir','Reverse')
% drawnow
% pause(.2)


%% comparison by range rating (would be good for future measurements)





