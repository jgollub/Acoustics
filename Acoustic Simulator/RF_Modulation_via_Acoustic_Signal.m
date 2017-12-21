%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Simulate Acoustic modulation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear
close all; clear; clc;
% Path
addpath('D:\Users\ho25\Desktop\Virtualizer\mexGpuFunctions\x64\Release');

%% Create monostatic measurement

% Configuration
% nProbes = 100;
nFreqs  = 1;
f   = linspace(20e9, 20e9, nFreqs);

tmp_probe = create_panel('type','probe_default','fsweep',f);
tmp_probe.u  =[0; 1; 0];
tmp_probe.v  =[0; 0; 1];

% probe array
probe_y_range=0.4;
probe_z_range=0.4;
sample_rate=.005;

yy = -probe_y_range/2:sample_rate:probe_y_range/2;
zz = -probe_z_range/2:sample_rate:probe_z_range/2;
xx = 0;

probes = cell(numel(zz),numel(yy),1);

in=0;
for i=1:numel(zz)
    for j=1:numel(yy)
        probes{i,j} = panel_offset(tmp_probe, 0, yy(j), zz(i));
%         probes{in} = panel_rotate(probes{in}, -(yy(j)/R), 0, 0, [R, 0, 0]);
    end
end

num_probes=numel(probes(:));
fprintf('The number of probes is %e \n',num_probes);
layout_fig=figure();
cla
 panel_plot(layout_fig,probes{1:end});

%% Setup model: wall in front of target object

target_res   = 0.0025;

% Load target
targetOffs=[];
targetOffs.x = 1;
targetOffs.y = 0;
targetOffs.z = 0;
targetOffs.alpha =  0*pi/180;
targetOffs.beta  =  0*pi/180;
targetOffs.gamma =  0*pi/180;

object  = objectCreator('Circle',1,2.5,15);
object  = objectMover(object,targetOffs); % Move
target  = ZBuffer(object,target_res); % Zbuffer

% Load wall
targetOffs=[];
targetOffs.x = .2;
targetOffs.y = 0;
targetOffs.z = 0;
targetOffs.alpha =  0*pi/180;
targetOffs.beta  =  0*pi/180;
targetOffs.gamma =  0*pi/180;

object  = objectCreator('Plate',1, 2.5, 2.5,40,40);
object  = objectMover(object,targetOffs); % Move
wall    = ZBuffer(object,target_res); % Zbuffer


% target0 = ZBuffer(object,2*target_res);

figure(layout_fig);
% hold on; scatter(target.locs(:,2),target.locs(:,3),5,'r'); 
hold on; scatter3(target.locs(:,1),target.locs(:,2),target.locs(:,3),5,'r'); 
hold on; scatter3(wall.locs(:,1),wall.locs(:,2),wall.locs(:,3),5,'r'); 

xlabel('x');
ylabel('y');
zlabel('z');
axis equal

%% Generate g: green's function response
opts.max_size = 2e5;
% num_voxels=length(kinect.locs);
% f = zeros(num_voxels,numel(probes));


mark=clock;
mark2=mark;
percent_check=1; last=0;
measurements_target=zeros(size(probes));
measurements_wall=zeros(size(probes));
for ii=1:size(probes,1)
    for jj =1:size(probes,2)
        
        measurements_target(ii,jj) = forward_model_SAR(probes{ii,jj},target.sigma,target.locs);
        measurements_wall(ii,jj)   = forward_model_SAR(probes{ii,jj},wall.sigma,wall.locs); 
  
        
        %tracking timing
        probNum=(floor(ii-1)*size(probes,2)+jj);
        if  (0 == mod(floor(probNum/num_probes*100), percent_check) & last~=floor(probNum/num_probes*100))
            last=floor(probNum/num_probes*100);
            
            time_step=etime(clock,mark2);
            time_step*(100-last)/percent_check;
            fprintf('finished %i%%. Time left is estimated to be %.1f Mins \n\n', last, ...
                (time_step*(100-last)/percent_check)/60);
            mark2=clock;
        end
    end
end
    fprintf('the total calculation time is, %4.3d seconds \n', etime(clock,mark))

%sum rows to get solution
% fsar=sum(f,2);
%% combine measurement



acoustic.f0=200;
acoustic.f1=200;
acoustic.t = 0:0.1E-3:10e-3; %milliseconds
acoustic.dt=t(2)-t(1);
acoustic.T=t(end)-t(1);

f_modulation=@(t) eps*2*pi*acoustic.f0*t;

acoustic.signal=real(exp(1i*f_modulation(acoustic.t))).';

fig_totalRF=figure();
fig_extracted=figure();
fig_ExtractedRF=figure();

measurements=zeros([size(measurements_target) length(t)]);
% 
% for eps=.5:.1:1.5;

% measurements=repmat(measurements_target,1,1,length(acoustic.t))...
%              +bsxfun(@times,measurements_wall,0.01*permute(exp(1.0i*f_modulation(t)),[1,3,2]));

measurements=repmat(measurements_target,1,1,length(acoustic.t))...
             +bsxfun(@times,measurements_wall,permute(acoustic.signal,[2,3,1]));

         
         
         %% apply Dan's extraction method to extract the target signal
Aprime=struct();
[Aprime.X,Aprime.Y]=meshgrid(yy,zz);
Aprime.f=f(1);
Aprime.measurements=removeModulatedRFSignal(measurements,f_modulation,acoustic.t,yy,zz);

%% reconstruct
singleFreqData=struct();
singleFreqData=Aprime;
nn=1;
%at measurement lane
singleFreqData.f=f(nn);


% plotSelectFields(fig_totalRF, 'up','fields',singleFreqData,nn);

%calc kspace for
singleFreqData = kComponents(singleFreqData,'measurements',nn);

% plotSelectFields(fig_totalRF, 'up','kspace',singleFreqData,nn);

c=3e8;
k0=2*pi*singleFreqData.f/c;
kzs=real(sqrt((2*k0).^2-singleFreqData.kx.^2-singleFreqData.ky.^2));

% zd=.1:.005:1.5;
% for ii=1:length(zd)
% Exy=(ifft2(fftshift(Ekxky.*exp(1j*kzs*zd(ii))))); %wave is traveling into -z henze positive phase to reverse
% 
% imagesc(X(1,:),Y(:,1),abs(Exy));
% axis equal; axis tight; axis xy; set(gca,'XDir','Reverse')
% title(['total fields, xd= ', num2str(zd(ii),'%.3f')])
% drawnow
% pause(.02)
% end
% figure(12)

figure(3);
d=1;
singleFreqData.ExyPropagated=(ifft2(fftshift(singleFreqData.Ekxky_measurements.*exp(1j*kzs*d)))); %wave is traveling into -z henze positive phase to reverse

imagesc(singleFreqData.Xupsampled(1,:),singleFreqData.Yupsampled(:,1),angle(singleFreqData.ExyPropagated))
axis equal; axis tight; axis xy;
title(['xd= ', num2str(d,'%.3f')])
set(gca,'XDir','Reverse')
drawnow
pause(.2)

%% Extract Target fields


% plotSelectFields(fig_ExtractedRF, 'down','fields',Aprime,nn);

Aprime=kComponents(Aprime,'measurements');

k0=2*pi*Aprime.f(nn)/c;
kzs=real(sqrt((2*k0).^2-Aprime.kx.^2-Aprime.ky.^2));

% plotSelectFields(fig_ExtractedRF, 'down','kspace', Aprime, nn);

% zd=.1:.005:1.5;
% for ii=1:length(zd)
% Exy=(ifft2(fftshift(Aprime.*exp(1j*kzs*zd(ii))))); %wave is traveling into -z henze positive phase to reverse
% 
% imagesc(X(1,:),Y(:,1),abs(Aprime));
% axis equal; axis tight; axis xy; set(gca,'XDir','Reverse')
% title(['Extracted fields, xd= ', num2str(zd(ii),'%.3f')])
% drawnow
% pause(.02)
% end

figure(4);
d=1;
Aprime.PropagatedAxy=(ifft2(fftshift(Aprime.Ekxky_measurements.*exp(1j*kzs*d)))); %wave is traveling into -z henze positive phase to reverse

imagesc(Aprime.Xupsampled(1,:),...
        Aprime.Yupsampled(:,1),abs(Aprime.PropagatedAxy))
axis equal; axis tight; axis xy;
title(['xd= ', num2str(d,'%.3f')])
set(gca,'XDir','Reverse')
drawnow
pause(.2)

end

%% add modulation

% subtract modulation using Dan's approach

% plot using Dan's approach
