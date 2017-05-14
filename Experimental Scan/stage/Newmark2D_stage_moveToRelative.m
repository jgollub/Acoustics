function [xUserPos, yUserPos] = Newmark2D_stage_moveToRelative(objg,speedmms,defZeroInXsteps,defZeroInYsteps,moveXmm,moveYmm)
steptomm=5000;
stepperspeed=round(speedmms*steptomm);
stepperaccel=round(10*stepperspeed);

CRLF =[char(13), char(10)];
%load home program
stagemove=['#MOVEREL', CRLF,...
'AC ',num2str(stepperaccel),',',num2str(stepperaccel), CRLF,...
'DC ',num2str(stepperaccel),',',num2str(stepperaccel), CRLF,...
'SP ',num2str(stepperspeed),',',num2str(stepperspeed), CRLF,...
'NOTE move relative position', CRLF,...
'PR ',num2str(-moveXmm*steptomm),',',num2str(-moveYmm*steptomm), CRLF,... % note negative sign used consistency with NFS coord sys
'BG', CRLF,...
'NOTE wait for end of move', CRLF,...
'AM', CRLF,...
'posx=_RPA', CRLF,...
'posy=_RPB', CRLF,...
'EN', CRLF];

objg.programDownload(stagemove);
response=objg.command(['XQ']);

testdone=1;
while (testdone>-1) || isnan(testdone)
    pause(.01);
    testdone=str2double(strtok(objg.command(['MG _XQ'])));
end
%xrange=strtok(objg.command(['MG rangex',CRLF]))
xposSteps=strtok(objg.command(['MG posx']));
xposSteps=str2double(xposSteps);

yposSteps=strtok(objg.command(['MG posy']));
yposSteps=str2double(yposSteps);

xUserPos=-(xposSteps-defZeroInXsteps)/steptomm; %note neg sign to match NFS coord sys
yUserPos=-(yposSteps-defZeroInYsteps)/steptomm;

end

