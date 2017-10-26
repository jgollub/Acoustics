function [ output_args ] = plotSelectFields(fHandle,location, type, measurementStruct,nn)

X            = measurementStruct.X;
Y            = measurementStruct.Y;
measurements = measurementStruct.measurements(:,:,nn);
f            = measurementStruct.f(nn);

if isfield(measurementStruct,'Ekxky')
    Ekxky        = measurementStruct.Ekxky;
    kx           = measurementStruct.kx;
    ky           = measurementStruct.ky;
end

ss = get(0,'screensize'); %The screen size
width = ss(3);
height = ss(4);
vert = 300; %300 vertical pixels
horz = 600; %600 horizontal pixels

figure(fHandle);
if strcmp('up',location)
    %This will place the figure in the top-right corner
    set(fHandle,'Position',[5/8*width-horz/2, height-vert, horz, vert]);
elseif strcmp('down',location)
    set(fHandle,'Position',[5/8*width-horz/2, 40, horz, vert]);
else
    disp('No position set');
end

switch type
    case 'fields'
        subplot(2,3,1)
        hold off;
        imagesc(X(1,:),Y(:,1),abs(measurements)); axis equal; axis tight; axis xy;
        set(gca,'XDir','Reverse')
        title('Abs(p)')
        
        subplot(2,3,2)
        hold off;
        imagesc(X(1,:),Y(:,1),real(measurements)); axis equal; axis tight; axis xy;
        set(gca,'XDir','Reverse')
        title('Real(p)')
        
        subplot(2,3,3)
        hold off;
        imagesc(X(1,:),Y(:,1),imag(measurements)); axis equal; axis tight;axis xy;
        set(gca,'XDir','Reverse')
        title('Imag(p)')
        hold on;
        
    case 'kspace'
        subplot(2,3,4)
        hold off;
        imagesc(kx(1,:),ky(:,1), abs(Ekxky)); axis tight; axis equal; axis xy;
        set(gca,'XDir','Reverse')
        title(['k-space af f= ',num2str(f)]);
        title('K-Space Abs(P)')
        
        subplot(2,3,5)
        hold off;
        imagesc(kx(1,:),ky(:,1), real(Ekxky)); axis tight; axis equal; axis xy;
        set(gca,'XDir','Reverse')
        title('K-Space Real(P)')
        
        subplot(2,3,6)
        hold off;
        imagesc(kx(1,:),ky(:,1), imag(Ekxky)); axis tight; axis equal; axis xy;
        set(gca,'XDir','Reverse')
        title('K-Space Imag(P)')
end   


