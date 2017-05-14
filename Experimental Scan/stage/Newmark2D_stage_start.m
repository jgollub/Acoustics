function objg = Newmark2D_stage_start()

%set return character
CRLF=[char(13), char(10)];

%create controller object for stage
objg = actxserver('galil');%set the variable g to the GalilTools COM wrapper
response = objg.libraryVersion;%Retrieve the GalilTools library versions
disp(response);%display GalilTools library version

objg.address = 'COM3 19200';%Open connections dialog box
response = objg.command(strcat(char(18), char(22)));%Send ^R^V to query controller model number
disp(strcat('Connected to: ', response));%print response

response = objg.command('MG_BN');%Send MG_BN command to query controller for serial number
disp(strcat('Serial Number: ', response));%print response

end

