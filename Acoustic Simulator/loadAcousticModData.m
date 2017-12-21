function [Struct] = loadAcousticModData(path)
Struct = struct();
load(path,'data');
Struct = data;

Struct.X=Struct.X/1000; 
Struct.dx=Struct.X(1,2)-Struct.X(1,1); %make sure dx is positive
Struct.Lx=abs(Struct.X(1,end)-Struct.X(1,1));

Struct.Y=Struct.Y/1000; 
Struct.dy=Struct.Y(2,1)-Struct.Y(1,1); %make sure d is positive
Struct.Ly=abs(Struct.Y(end,1)-Struct.Y(1,1)); 

% Struct.k0=2*pi*Struct.f(nn)/c;

Struct.dt              = Struct.sweepTime/(Struct.numPoints-1);
Struct.CWSamplingFreqs = linspace(0,1/Struct.dt,Struct.numPoints);
Struct.tMeas  = linspace(0,Struct.sweepTime,Struct.numPoints);

