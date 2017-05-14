function [xposSteps,yposSteps] = Newmark2D_stage_getStepPosition(objg)
steptomm=5000;
xposSteps=str2double(strtok(objg.command('MG _RPA')));
yposSteps=str2double(strtok(objg.command('MG _RPB')));
end