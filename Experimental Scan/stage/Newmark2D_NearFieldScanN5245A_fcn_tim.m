function [X,Y,f,measurements] = Newmark2D_NearFieldScanN5245A_fcn_tim...
    (objg,speedmms,defZeroInXsteps,defZeroInYsteps,...
    xmin,xmax,ymin,ymax,dstep,frequencysamples,...
    fstart,fstop,switches,IFbandwidth,calfile,~)


if switches ~= 256 && switches~=0 && switches~=2
    error('not the correct number of switches/modes')
end

%IF bandwidth is set to
IFbandwidth
%cal file is set to
calfile
%select measurment
PortsToMeas='MeasS31';


%% initialize instruments
delete(instrfind) %delete any existing instrument objects 

%%%%%%%%%%%%% INSERT ARDUINO MECHANISM %%%%%%%%%%%%%%%%%%
if switches>0
    s = startSerial('COM4');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vobj_vna = agilent_N5245A_NA_startVISAObject_TCPIP;     % open vna communications

[buffersize, f] = Initialize_N5245A_FastwithThruCal(...
    vobj_vna, frequencysamples, fstart, fstop, IFbandwidth, calfile); % setup VNA scan

%% setup scan 
[X, Y] = meshgrid(xmin:dstep:xmax,ymin:dstep:ymax);
if switches>0
    measurements = zeros(size(Y,1),size(X,2),frequencysamples,switches); %changed size(X,1) to size(X,2)
else
    measurements = zeros(size(Y,1),size(X,2),frequencysamples); %changed size(X,1) to size(X,2)
end
stops = size(Y,1)*size(X,2);

%% begin scan!
tic
stopscomp = 0;

for yn=1:size(Y,1)
    direction = 2*mod(yn,2)-1;
    if direction>0
        xindex = 1:size(X,2);
    else
        xindex = size(X,2):-1:1;
    end
    for xn=xindex
        x = X(yn,xn); %no flipdim here--we are interested flipped position anyway
        y = Y(yn,xn); %to account for meshgrid making array that starts at neg pos and goes pos
        Newmark2D_stage_moveToAbsolute(objg,speedmms,defZeroInXsteps,defZeroInYsteps,x,y); %recommended speed is 25 mm/sec
                
        %%%%%%%%%%%%% INSERT ARDUINO MECHANISM %%%%%%%%%%%%%%%%%% 
        if switches==2
            for jj=1:2
                aps=[0 255];
                %shiftOut(ard,8,2,aps(jj))
                dataString = strcat(num2str(1),',',num2str(aps(jj)));
                fprintf(s,dataString);
                for qq=1:20000000; end
                measurements(yn,xn,:,jj) = Read_N5245A(vobj_vna,buffersize,PortsToMeas);
            end
        elseif switches>2
            for jj=1:switches
                aps=0:255;
                %shiftOut(ard,8,2,jj)
                dataString = strcat(num2str(1),',',num2str(aps(jj)));
                fprintf(s,dataString);
                for qq=1:20000000; end
                measurements(yn,xn,:,jj) = Read_N5245A(vobj_vna,buffersize,PortsToMeas);
            end
        else
            measurements(yn,xn,:) = Read_N5245A(vobj_vna,buffersize,PortsToMeas);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        stopscomp = stopscomp+1;
        timere = (stops-stopscomp)*toc;
        disp(['Est. time remaining: ' num2str(timere/60) ' minutes'])

        tic  
    end
    
    figure(2)
    if switches>0
    imagesc(abs(measurements(:,:,round(length(measurements(1,1,:,1))/2),1)))
    else
    imagesc(abs(measurements(:,:,round(length(measurements(1,1,:,1))/2))))
    end
end

for hide=1; %send you a friendly email
    
data.X = X;
data.Y = Y;
data.f = f;
data.measurements = measurements;
save(['C:\Users\Lab\Desktop\3-D Mapper Code\results\NFS-' date],'data')

dataBytes=whos('data');
Dbytes=dataBytes.bytes;

myaddress = 'VNA2YOU@gmail.com';
mypassword = '6plasmons';

setpref('Internet','E_mail',myaddress);
setpref('Internet','SMTP_Server','smtp.gmail.com');
setpref('Internet','SMTP_Username',myaddress);
setpref('Internet','SMTP_Password',mypassword);

props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
    'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
    
message1=['Hey buddy!' 10 10 'I was wondering...' 10 10 'Do you know what the definition of a polar bear is?' 10 10 'Its a rectangular bear after a coordinate transform!' 10 10 10 'Also, I wanted to let you know that I am all done and you can get your data whenever you are ready! Look forward to seeing you!' 10 10 'Your friend,' 10 'VNA + NFS' 10 10 'P.S. - I am attaching the data in case you are too lazy to walk over here...'];
message2=['Hey buddy!' 10 10 'I was wondering...' 10 10 'Do you know what the definition of a polar bear is?' 10 10 'Its a rectangular bear after a coordinate transform!' 10 10 10 'Also, I wanted to let you know that I am all done and you can get your data whenever you are ready! Look forward to seeing you!' 10 10 'Your friend,' 10 'VNA + NFS' 10 10 'P.S. - Sorry, your file was too big or I would have sent it to you...'];

for recipients={'sleasmant@gmail.com'}
if Dbytes<12E6;
    sendmail(recipients,'Your near-field scan is ready!',message1,['C:\Users\Lab\Desktop\3-D Mapper Code\results\NFS-' date '.mat'])
else
    sendmail(recipients,'Your near-field scan is ready!',message2)
end
end

end

%% clean up communications
fclose(s);
fprintf(vobj_vna,'DISP:ENABLE on');
agilent_N5245A_NA_stopVISAObject(vobj_vna);
delete(instrfind);

