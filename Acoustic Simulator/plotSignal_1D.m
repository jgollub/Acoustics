%check frequency domain of signal

measureStruct=WithAcousticMod;
nx=50;
ny=50;

%plot raw fields
figure(1)
subplot(5,2,1)
plot(db(squeeze(squeeze(measureStruct.measurements(ny,nx,:)))),'k');
axis tight;

subplot(5,2,2)
plot(db(squeeze(squeeze(measureStruct.measurements(ny,nx,:)))),'r');
axis tight;

subplot(5,2,3:6)
diff=db(squeeze(squeeze(measureStruct.measurements(ny,nx,:)...
                     -WithoutAcousticMod.measurements(ny,nx,:))));
plot(diff)
axis tight;

%frequency sampling (for VNA measurements)
acousticFFTfreqs=0:1/T:1/dt;
measurement_fft=fft(squeeze(squeeze(measureStruct.measurements(ny,nx,:))))

subplot(5,2,7:10)
plot(acousticFFTfreqs,db(measurement_fft))
axis tight;
title('fft signal')
xlabel('frequency (Hz)')


%% plot 200 Hz magnitude for each position
measurement_fft=zeros([size(measureStruct.X),length(measureStruct.f)]);

%pick out acoustic signal
acousticFreqs=200;
deltaFreq=5;


for nx = 1:length(measureStruct.X(1,:))
    for ny = 1:length(measureStruct.Y(:,1))
        
        measurement_fft(ny,nx,:) = squeeze(squeeze(fft(measureStruct.measurements(nx,ny,:))));
        
    end
end

[minVal, minIndx]=find(deltaFreq > abs(acousticFreqs-acousticFFTfreqs));

% SumAtAcousticFields=sum(measurement_fft(:,:,minIndx),3) -mean(mean(mean(measurement_fft)));
SumAtAcousticFields=sum(measurement_fft(:,:,minIndx),3);


figure()
imagesc(X(1,:),Y(:,1),abs(SumAtAcousticFields))
axis equal; axis tight; axis xy;
title(['xd= ', num2str(d,'%.3f')])
set(gca,'XDir','Reverse')
drawnow
pause(.2)





