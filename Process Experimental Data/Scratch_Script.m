t=0:1:1000;
modulated_signal = 10*sin(2*pi*(1/100)*t)+10
figure(100); clf;
subplot(2,1,1)
plot(t,modulated_signal)
title('signal')
subplot(2,1,2)

plot(xcorr(modulated_signal,modulated_signal))
title('cross-correlation')