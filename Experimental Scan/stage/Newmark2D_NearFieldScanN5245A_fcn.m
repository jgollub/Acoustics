

function [X,Y,f,measurements] = NearFieldScanN5245A_fcn(objg,speedmms,defZeroInXsteps,defZeroInYsteps,...
    xmin,xmax,ymin,ymax,dstep,frequencysamples,fstart,fstop,switches,IFbandwidth,calfile,measureSparam)

%IF bandwidth is set to
IFbandwidth

%cal file is set to
calfile

%select measurment
if measureSparam==1
    PortsToMeas='MeasS31'
elseif measureSparam==2
    PortsToMeas='MeasS32'
elseif measureSparam==3
    PortsToMeas='MeasS21'
end

%% initialize instruments
delete(instrfind) %delete any existing instrament objects 
if switches
    vobj_switch = agilent_11713C_switchdriver_startVISAObject; %open switch communications
end
vobj_vna = agilent_N5245A_NA_startVISAObject_TCPIP;     % open vna communications

%IFbandwidth=1000;
[buffersize, f] = Initialize_N5245A_FastwithThruCal(vobj_vna, frequencysamples, fstart, fstop, IFbandwidth, calfile); % setup VNA scan
%[buffersize, f] = Initialize_N5245A_vers3(vobj_vna, frequencysamples, fstart, fstop, IFbandwidth, calfile); % setup VNA scan
%[buffersize, f] = Initialize_8364B(vobj_vna, frequencysamples, fstart, fstop); % setup VNA scan

%% setup scan 
[X, Y] = meshgrid(xmin:dstep:xmax,ymin:dstep:ymax);
%[X, Y] = meshgrid(xmin:dstep:xmax,ymax:-dstep:ymin);

switchbank=1;
if switches>6
    switchbank=2;
end
if switches
    agilent_11713C_switchdriver_closeAll(vobj_switch) %pin 16 [open all] (LED button number 9) must be on for -024 option switches 
    %agilent_11713C_switchdriver_openAll(vobj_switch)
    for sb=1:switchbank
        agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,sb,[1:(switches/switchbank)])
     %  agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,sb,[1:(switches/2)]) %sloppy... accounting for two banks of
    end
    measurements = zeros(size(Y,1),size(X,2),frequencysamples,switches); %changed size(X,1) to size(X,2)
else
    measurements = zeros(size(Y,1),size(X,2),frequencysamples); %changed size(X,1) to size(X,2)
end

stops = size(Y,1)*size(X,2);

%% begin scan!
tic
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
        x = X(yn,xn); %no flipdim here--we are interested flipped position anyway
        %y = Y(yn,xn); 
        y = Y(yn,xn); %to account for meshgrid making array that starts at neg pos and goes pos
        
        Newmark2D_stage_moveToAbsolute(objg,speedmms,defZeroInXsteps,defZeroInYsteps,x,y); %recommended speed is 25 mm/sec
        
        %are we to use the switch?
        if switches %if yes...
            agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,1,[1:switches/2])
            for sb=1:switchbank
                for sn=1:(switches/switchbank)
                    if sb==1
                        dummymeasure=PortsToMeas;
                    elseif sb==2
                        dummymeasure='MeasS31';
                    end
                    agilent_11713C_switchdriver_openChannelNumbers(vobj_switch,sb,sn)
                    %collect data from VNA
                    measurements(yn,xn,:,(sb-1)*(switches/switchbank)+sn) = Read_N5245A(vobj_vna,buffersize,dummymeasure);
                    %                 measurements(yn,xn,:,sn) = Read_8364B(vobj_vna,buffersize);
                    agilent_11713C_switchdriver_closeChannelNumbers(vobj_switch,sb,sn)
                end
            end
        else %if not
            %collect data from VNA
            %!!!!!!!!!!!!!!! need to fix GUI to select correct measurement
            measurements(yn,xn,:) = Read_N5245A(vobj_vna,buffersize,'MeasS31');
%           measurements(yn,xn,:) = Read_8364B(vobj_vna,buffersize);
        end        
        stopscomp = stopscomp+1;
        timere = (stops-stopscomp)*toc/3600;
        disp(['Est. time remaining: ' num2str(timere) 'hours'])
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
    end
    if switches>6
        figure(2)
        imagesc(abs(measurements(:,:,round(length(measurements(1,1,:,1))/2),1)))
        figure(3)
        imagesc(abs(measurements(:,:,round(length(measurements(1,1,:,1))/2),7)))
    elseif (switches<=6) && (switches>0)
        figure(2)
        imagesc(abs(measurements(:,:,round(length(measurements(1,1,:,1))/2),1)))
    else
        figure(2)
        imagesc(abs(measurements(:,:,round(length(measurements(1,1,:))/2))))
    end
end

%% clean up communications
fprintf(vobj_vna,'DISP:ENABLE on');
if switches
    agilent_11713C_switchdriver_stopVISAObject(vobj_switch);
end
%agilent_E8364B_NA_stopVISAObject(vobj_vna);
agilent_N5245A_NA_stopVISAObject(vobj_vna);
delete(instrfind);