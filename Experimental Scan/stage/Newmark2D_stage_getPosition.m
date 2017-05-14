function [xUserPos,yUserPos] = Newmark2D_stage_getPosition(objg,defZeroInXsteps,defZeroInYsteps)
steptomm=5000;
xposSteps=str2double(strtok(objg.command('MG _RPA')));
yposSteps=str2double(strtok(objg.command('MG _RPB')));

xUserPos=-(xposSteps-defZeroInXsteps)/steptomm; %note negative sign to match NFS coord Sys
yUserPos=-(yposSteps-defZeroInYsteps)/steptomm;
end

