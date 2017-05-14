function [] = Newmark2D_stage_StopStage(objg)

%set return character
CRLF =[char(13), char(10)];

% stop stage immediately
response=objg.command('ST');

end

