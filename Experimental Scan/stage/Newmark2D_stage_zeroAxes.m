function [defZeroInXsteps, defZeroInYsteps]=Newmark2D_stage_zeroAxes(objg)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

defZeroInXsteps=str2double(strtok(objg.command('MG _RPA')));
defZeroInYsteps=str2double(strtok(objg.command('MG _RPB')));

end

