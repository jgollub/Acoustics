function [] = plotPropgatedFields(fHandle,location,fieldStruct)

ss = get(0,'screensize'); %The screen size
width = ss(3);
height = ss(4);
vert = 300; %300 vertical pixels
horz = 400; %600 horizontal pixels

figure(fHandle);
if strcmp('up',location)
    %This will place the figure in the top-right corner
    set(fHandle,'Position',[10.75/12*width-horz/2, height-vert, horz, vert]);
elseif strcmp('down',location)
    set(fHandle,'Position',[10.75/12*width-horz/2, 40, horz, vert]);
else
    disp('No position set');
end

imagesc(fieldStruct.X(1,:),fieldStruct.Y(:,1),abs(fieldStruct.PropagatedExy))
axis equal; axis tight; axis xy;
title(['xd= ', num2str(fieldStruct.d,'%.3f')])
set(gca,'XDir','Reverse')
drawnow
