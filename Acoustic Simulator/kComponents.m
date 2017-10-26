function [measurementStruct] = kComponents(measurementStruct,nn)

X            = measurementStruct.X;
Y            = measurementStruct.Y;

if length(size(measurementStruct.measurements)) == 3;
    measurements = measurementStruct.measurements(:,:,nn);
else
    measurements = measurementStruct.measurements;
end

dx=X(1,2)-X(1,1);
dy=Y(2,1)-Y(1,1);

[xnum,ynum]=size(X);

pad=2^nextpow2(max(size(measurements)));
[kx,ky]=meshgrid(linspace(-2*pi/(2*dx),2*pi/(2*dx),pad),linspace(-2*pi/(2*dy),2*pi/(2*dy),pad));
dkx=kx(1,2)-kx(1,1);
dky=ky(2,1)-ky(1,1);

shiftx=dx*(pad/2-ceil(xnum/2));
shifty=dy*(pad/2-ceil(ynum/2));

Ekxky=circshift(fft2(measurements, pad, pad),[pad/2,pad/2]);
Ekxky=Ekxky.*exp(-1.0j*kx*(shiftx)).*exp(-1.0j*ky*(shifty));

measurementStruct.Ekxky      = Ekxky;
measurementStruct.kx         = kx;
measurementStruct.ky         = ky;
[measurementStruct.Xupsampled, ...
 measurementStruct.Yupsampled] = meshgrid(linspace(-2*pi/(2*dkx),2*pi/(2*dkx), pad),...
                                          linspace(-2*pi/(2*dky),2*pi/(2*dky), pad));
end

