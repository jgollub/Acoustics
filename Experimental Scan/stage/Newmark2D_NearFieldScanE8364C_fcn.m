function [data] = Newmark2D_NearFieldScanE8364C_fcn(objg,speedmms,defZeroInXsteps,defZeroInYsteps,...
    xmin,xmax,ymin,ymax,dstep,NumFreqs,fstart,fstop,NumApMasks,IFBW,calfile,~)

%% enter savename, aperture masks, and display the settings for user
power=3;
savename='PEC_camo_24_3';
sParMeas='S21';

fprintf('THIS IS FROM TIMAGER DIRECTORY \n')

% ADD HARDWARE TO THE PATH

%NOTE!! 
%NUMBER OF SHIFT REGISTERS IS HARDCODED AS 14 
%ARD SKETCH MUST CHANGE IF THIS IS CHANGED
%MAKE SURE THE ARDUINO IS SET UP ON THE CORRECT COM PORT

if 0
load('A_masks.mat');
ApMasks1=1-A_masks;

m1=576;
end

%take a binary matrix with m1 # 
%of modes and make it a decimal version
% ApMasks=reshape(ApMasks1',8,m1*14)';
% ApMasks=bi2de(ApMasks,'left-msb');
% ApMasks=reshape(ApMasks,14,m1)';
% clearvars custModes eyeModes m1 sh t1l t1r

if NumApMasks==0
    fprintf('\nNo masks are being used\n')
    ApMasks=[];
elseif size(ApMasks,1)~=NumApMasks
    error('The number of GUI switches does not match the data hard-coded into the _fcn')
else
end

s=sprintf(['\nIFBW = %g\nCal = ',calfile,...
    '\nnumber of masks = %i \npower = %i dBm\nMeasuring ',sParMeas,' by hard-coding\n'], IFBW,NumApMasks,power);
disp(s)

%% initialize instruments
delete(instrfind) %delete any existing instrument objects 

vobj_vna=InstrumentOpen('visa','agilent','E8364C');

[buffersize,f]=VNAInitiate(vobj_vna,NumFreqs,fstart,fstop,IFBW,power,sParMeas,calfile);

if NumApMasks>0
    serPort = startSerial('COM9'); %must be done after instrument open
end

%% setup scan 
[X, Y] = meshgrid(xmin:dstep:xmax,ymin:dstep:ymax);

if NumApMasks
    measurements = zeros(size(Y,1),size(X,2),NumFreqs,NumApMasks);
else
    measurements = zeros(size(Y,1),size(X,2),NumFreqs);
end

stops = size(Y,1)*size(X,2);

%% begin scan!
tic
stopscomp=0;

%take the measurement INSIDE this IF statement
if NumApMasks
    for yn=1:size(Y,1)
        direction = 2*mod(yn,2)-1;
        if direction>0
            xindex = 1:size(X,2);
        else
            xindex = size(X,2):-1:1;
        end
        for xn=xindex
            x = X(yn,xn);
            y = Y(yn,xn);
            Newmark2D_stage_moveToAbsolute(objg,speedmms,defZeroInXsteps,defZeroInYsteps,x,y); %recomm. speed 25 mm/sec
            
            for jj=1:NumApMasks
                ApTemp=ApMasks(jj,:);
                str_x = strcat(num2str(14),',',num2str(ApTemp));
                fprintf(serPort,str_x);
                measurements(yn,xn,:,jj) = VNARead(vobj_vna,buffersize,sParMeas);
            end
            
            stopscomp = stopscomp+1;
            timere = (stops-stopscomp)*toc/3600;
            disp(['Est. time remaining: ' num2str(timere) ' hours'])
            tic
        end
        figure(2); imagesc(abs(measurements(:,:,round(length(measurements(1,1,:,1))/2),1)))
    end
else % ApMasks=0
    for yn=1:size(Y,1)
        direction = 2*mod(yn,2)-1;
        if direction>0
            xindex = 1:size(X,2);
        else
            xindex = size(X,2):-1:1;
        end
        for xn=xindex
            x = X(yn,xn);
            y = Y(yn,xn);
            Newmark2D_stage_moveToAbsolute(objg,speedmms,defZeroInXsteps,defZeroInYsteps,x,y); %recomm. speed 25 mm/sec
            
            measurements(yn,xn,:) = VNARead(vobj_vna,buffersize,sParMeas);
            
            stopscomp = stopscomp+1;
            timere = (stops-stopscomp)*toc/3600;
            disp(['Est. time remaining: ' num2str(timere) ' hours'])
            tic
        end
        figure(2); imagesc(abs(measurements(:,:,round(length(measurements(1,1,:))/2))))
    end
end

data.X = X;
data.Y = Y;
data.f = f;
data.masks=ApMasks;
data.measurements = measurements;
a=clock;
save(['C:\Users\Floor 1 Imager\Desktop\NFS data\camo\',savename,'_',date,'_',num2str(a(4)),'_',num2str(a(5))],'data')

%% put it in an email
if 0
dataBytes=whos('data'); Dbytes=dataBytes.bytes;

myaddress = 'VNA2YOU@gmail.com'; mypassword = '6plasmons';

setpref('Internet','E_mail',myaddress); setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress); setpref('Internet','SMTP_Password',mypassword);

props = java.lang.System.getProperties; props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
    'javax.net.ssl.SSLSocketFactory'); props.setProperty('mail.smtp.socketFactory.port','465');

message=['Dear NFS user,' 10 10 ...
         'The deed is done.' 10 10 ...
         '-Your friendly neighborhood NFS'];
        
 for recipients={'mboyarsky12@gmail.com','laura.pulido.mancera@gmail.com','sleasmant@gmail.com'}
    if Dbytes<12E6;
        sendmail(recipients,'NFS is back in action!',message,...
        ['C:\Users\Floor 1 Imager\Desktop\NFS data\beamsteering\',savename,'_',date,'_',num2str(a(4)),'_',num2str(a(5)),'.mat']) %att. file
    else
        sendmail(recipients,'NFS is back in action!',message)
    end
 end
end

%% clean up communications
fprintf(vobj_vna,'DISP:ENABLE on');
InstrumentClose();