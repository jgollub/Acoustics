function [xrange,yrange] = Newmark2D_stage_RangeOfStage(objg,speedmms)
steptomm=5000;
stepperspeed=round(speedmms*steptomm);
stepperaccel=round(10*stepperspeed);
%set return character
CRLF =[char(13), char(10)];
%load home program
stagerangetest=['#RANGEST', CRLF,...
'AC ',num2str(stepperaccel),',',num2str(stepperaccel), CRLF,...
'DC ',num2str(stepperaccel),',',num2str(stepperaccel), CRLF,...
'JG ',num2str(-stepperspeed),',',num2str(-stepperspeed), CRLF,...
'NOTE jog until you hit limits', CRLF,...
'BG', CRLF,...
'NOTE wait for end of move', CRLF,...
'AM', CRLF,...
'NOTE take 1mm step back from edge', CRLF,...
'PR 5000,5000', CRLF,...
'BG', CRLF,...
'AM', CRLF,...
'NOTE capture min,max step val', CRLF,...
'minx=_RPA', CRLF,...
'miny=_RPB', CRLF,...
'NOTE jog to other side and record pos', CRLF,...
'JG ',num2str(stepperspeed),',',num2str(stepperspeed), CRLF,...
'BG', CRLF,... 
'NOTE jog to the other side', CRLF,...
'AM', CRLF,...
'maxx=_RPA', CRLF,...
'maxy=_RPB', CRLF,...
'NOTE find range in steps', CRLF,...
'rangex=maxx-minx', CRLF,...
'rangey=maxy-miny', CRLF,...
'AM',CRLF,...
'NOTE fset back to standard stage spead acel/decel', CRLF,...
'SP 125000,125000', CRLF,...
'AC 1250000,1250000', CRLF,...
'DC 1250000,1250000', CRLF,...
'PR -5000,-5000' CRLF,...
'BG', CRLF,...
'AM',CRLF,...
'EN', CRLF];

objg.programDownload(stagerangetest);
response=objg.command(['XQ']);

%wait for program to finish

testdone=1;
while (testdone>-1) || isnan(testdone)
    pause(.01);
    testdone=str2double(strtok(objg.command(['MG _XQ'])));
end
%xrange=strtok(objg.command(['MG rangex',CRLF]))
xrange=strtok(objg.command(['MG rangex']));
xrange=str2double(xrange);

yrange=strtok(objg.command(['MG rangey']));
yrange=str2double(yrange);

xrange=(xrange/steptomm); %subtract 1mm to keep stage away from limit switch
yrange=(yrange/steptomm);

fileID = fopen('C:\Users\Lab\Documents\MATLAB\Newmark 2D Stage\RangeOfStage','w');
fprintf(fileID,'%7s %7s\r\n','X Range','Y Range');
fprintf(fileID,'%6.4f %6.4f\n',xrange,yrange);
fclose(fileID);
end
