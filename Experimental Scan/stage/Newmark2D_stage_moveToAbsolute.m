function [xUserPos,yUserPos] = Newmark2D_stage_moveToAbsolute(objg,speedmms,defZeroInXsteps,defZeroInYsteps,moveToUserXmm,moveToUserYmm)
steptomm=5000;
stepperspeed=round(speedmms*steptomm);
stepperaccel=round(10*stepperspeed);

CRLF =[char(13), char(10)];
%load home program
stagemove=['#MOVEABS', CRLF,...
'AC ',num2str(stepperaccel),',',num2str(stepperaccel), CRLF,...
'DC ',num2str(stepperaccel),',',num2str(stepperaccel), CRLF,...
'SP ',num2str(stepperspeed),',',num2str(stepperspeed), CRLF,...
'NOTE move to abs position', CRLF,...
'PA ',num2str(-moveToUserXmm*steptomm+defZeroInXsteps),',',num2str(-moveToUserYmm*steptomm+defZeroInYsteps), CRLF,... % note negative sign used consistency with NFS coord sys
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
xpos=strtok(objg.command(['MG posx']));
xpos=str2double(xpos);

ypos=strtok(objg.command(['MG posy']));
ypos=str2double(ypos);

xUserPos=-(xpos-defZeroInXsteps)/steptomm; %note neg sign to match NFS coord sys
yUserPos=-(ypos-defZeroInYsteps)/steptomm;

end

