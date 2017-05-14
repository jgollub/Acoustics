function [] = Newmark2D_stage_stop(objg)
%clear object
delete(objg);%delete g object and close connection with controller
disp('Controller connection stopped')
end

