
%% load data
load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\RFAcousticModulation\wAcoustic_withMetalPlate_wPost_NFS-29-Sep-2017.mat','data');
WithAcousticMod = data;
load('D:\Dropbox (Duke Electric & Comp)\MetaImager Data\Acoustic Data\RFAcousticModulation\PostInFront_wOutAcousticNFS-22-Sep-2017.mat','data');
WithoutAcousticMod = data;
clear data;

f=WithAcousticMod.f;
c=3e8;


X=WithAcousticMod.X/1000; 
dx=X(1,2)-X(1,1); %make sure dx is positive
Lx=abs(X(1,end)-X(1,1));
[ynum,xnum]=size(X); 

Y=WithAcousticMod.Y/1000; 
dy=Y(2,1)-Y(1,1); %make sure d is positive
Ly=abs(Y(end,1)-Y(1,1)); 

%%  compare modulated and unmodulated

fieldDifference=WithAcousticMod;
nn=1;
d=0.35;
for nn=1:length(WithAcousticMod.f);
fieldDifference.measurements=1*WithAcousticMod.measurements(:,:,nn)-1*WithoutAcousticMod.measurements(:,:,nn);
fieldDifference=kComponents(fieldDifference);

k0=2*pi*fieldDifference.f/c;
kzs=real(sqrt((2*k0).^2-fieldDifference.kx.^2-fieldDifference.ky.^2));
fieldDifference.PropagatedExy = ifft2(fftshift(fieldDifference.Ekxky.*exp(1j*kzs*d)));
fig_plotDiff = figure(10);

imagesc(X(1,:),Y(:,1),db(fieldDifference.PropagatedExy),[-60,-40]);
axis equal; axis tight; axis xy;
title(['xd= ', num2str(d,'%.3f')])
set(gca,'XDir','Reverse')
drawnow; disp(['step', num2str(nn)])

pause(.1)
end


%%

% propagate fields single frequency
nn=1;
singleFreqExp_wAcoust = WithAcousticMod.measurements(:,:,nn);

WithAcousticMod.X            = X;
WithAcousticMod.Y            = Y;

fig_ExpTotalRF=figure(20);
plotSelectFields(fig_ExpTotalRF, 'up','fields',WithAcousticMod,nn);

WithAcousticMod = kComponents(WithAcousticMod,nn);

plotSelectFields(fig_ExpTotalRF, 'up','kspace',WithAcousticMod,nn);

k0=2*pi*WithAcousticMod.f(nn)/c;
kzs=real(sqrt((2*k0).^2-WithAcousticMod.kx.^2-WithAcousticMod.ky.^2));

% figure(14);
% zd=.1:.005:.5;
% 
% for ii=1:length(zd)
% ExyPropagated=(ifft2(fftshift(singleFreqExp_wAcoust.Ekxky.*exp(1j*kzs*zd(ii))))); %wave is traveling into -z henze positive phase to reverse
% 
% imagesc(X(1,:),Y(:,1),abs(ExyPropagated));
% axis equal; axis tight; axis xy; set(gca,'XDir','Reverse')
% title(['xd= ', num2str(zd(ii),'%.3f')])
% drawnow
% pause(.05)
% end

fig_propagated=figure(21)
WithAcousticMod.d=.25;
WithAcousticMod.PropagatedExy=(ifft2(fftshift(WithAcousticMod.Ekxky.*exp(1j*kzs*WithAcousticMod.d)))); %wave is traveling into -z henze positive phase to reverse

plotPropgatedFields(fig_propagated,'up',WithAcousticMod)


% figure()
% 
% WithAcousticMod=kComponents(WithAcousticMod);
% 
% k0=2*pi*WithAcousticMod.f(nn)/c;
% kzs=real(sqrt((2*k0).^2-WithAcousticMod.kx.^2-WithAcousticMod.ky.^2));
% 
% plotSelectFields(fig_ExtractedRF, 'down' ,'kspace', WithAcousticMod,1);
% 
% fig_propagated=figure();
% WithAcousticMod.d=.35;
% WithAcousticMod.PropagatedExy=(ifft2(fftshift(WithAcousticMod.Ekxky.*exp(1j*kzs*WithAcousticMod.d)))); %wave is traveling into -z henze positive phase to reverse
% 
% plotPropgatedFields(fig_propagated,'down',WithAcousticMod)



%% Remove signal from wall using acoustic modulation
t  = linspace(0,WithAcousticMod.sweepTime,WithAcousticMod.numPoints);
dt = t(2)-t(1);
T  = WithAcousticMod.sweepTime;


if isfield(WithAcousticMod,'acousticSignal')
    f_a          = WithAcousticMod.acousticSignal.f0;
else
    f_a          = 0.2E3; %!!!!!!!!!!!!!!!!! HARD CODED 
end
f_modulation = @(t) 2*pi*f_a*t;

Aprime=struct();
Aprime.X = X;
Aprime.Y = Y;
Aprime.f=f;

Aprime.measurements=removeModulatedRFSignal(WithAcousticMod.measurements,f_modulation,t,WithAcousticMod.X,WithAcousticMod.Y);

%% Extract Target fields
fig_ExtractedRF=figure(22);
plotSelectFields(fig_ExtractedRF, 'down', 'fields',Aprime,1);

Aprime=kComponents(Aprime);

k0=2*pi*Aprime.f(nn)/c;
kzs=real(sqrt((2*k0).^2-Aprime.kx.^2-Aprime.ky.^2));

plotSelectFields(fig_ExtractedRF, 'down' ,'kspace', Aprime,1);

fig_propagated=figure();
Aprime.d=.35;
Aprime.PropagatedExy=(ifft2(fftshift(Aprime.Ekxky.*exp(1j*kzs*Aprime.d)))); %wave is traveling into -z henze positive phase to reverse

plotPropgatedFields(fig_propagated,'down',Aprime)

% imagesc(X(1,:),Y(:,1),abs(Aprime.PropagatedAxy))
% axis equal; axis tight; axis xy;
% title(['xd= ', num2str(d,'%.3f')])
% set(gca,'XDir','Reverse')
% drawnow
% pause(.2)






