function [] = plotPropgatedFields(fHandle,plotLoc,struct,structField_PropagatedExy,varargin)
p=inputParser;

%example inputs
%plotPropgatedFields(fHandle,location,Struct,PropagatedExy,varargin)
%
% 'fHandle'         figure handle
% 'plotLoc'         plot location
% 'Struct'          structure containing measured/simulated data
% 'PropagatedExy'   structure field that contains propagatedfields
%
% 'clim'            plot limits
% 'plotType'        plot form {db, real, imaginary, angle}

expectedPlotLoc={'up','down','sidebyside_left','sidebyside_right'};
expectedPlotType={@db, @real, @imaginary, @angle};


addRequired(p,'fHandle',@ishandle);
addRequired(p,'plotLoc',@(x) any(validatestring(x,expectedPlotLoc)));
addRequired(p,'struct',@isstruct);
addRequired(p,'structField_PropagatedExy',@(x) isfield(struct,x));

addParameter(p,'plotType', @db,@(x) isa(x,'function_handle'));
addParameter(p,'clim',NaN,@isnumeric);

parse(p,fHandle,plotLoc,struct,structField_PropagatedExy,varargin{:});

ss = get(0,'screensize'); %The screen size
width = ss(3);
height = ss(4);
vert = 300; %300 vertical pixels
horz = 400; %600 horizontal pixels

figure(p.Results.fHandle);
if strcmp('up',p.Results.plotLoc)
    %This will place the figure in the top-right corner
    set(p.Results.fHandle,'Position',[10.75/12*width-horz/2, height-vert, horz, vert]);
elseif strcmp('down',p.Results.plotLoc)
    set(p.Results.fHandle,'Position',[10.75/12*width-horz/2, 40, horz, vert]);
elseif strcmp('sidebyside_left',p.Results.plotLoc)
    set(p.Results.fHandle,'Position',[10.75/12*width-horz/2, 40+horz, 2*horz, vert]);
    subplot(1,2,1)
elseif strcmp('sidebyside_right',p.Results.plotLoc)
    set(p.Results.fHandle,'Position',[10.75/12*width-3*horz/2, 40+horz, 2*horz, vert]);
    subplot(1,2,2)
end

if any(p.Results.clim)
        imagesc(struct.X(1,:),struct.Y(:,1),...
        p.Results.plotType(struct.(structField_PropagatedExy)),...
        p.Results.clim);
    
else
    imagesc(struct.X(1,:),struct.Y(:,1),...
        p.Results.plotType(struct.(structField_PropagatedExy))...
        );
end

axis equal; axis tight; axis xy;
title(['xd= ', num2str(struct.d,'%.3f')])
set(gca,'XDir','Reverse')
drawnow
