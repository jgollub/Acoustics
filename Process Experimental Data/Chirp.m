
%plot input chirp
figure(100)
cla
plot(t,real(chirp_sig))
hold on
plot(t,real(real(signalOut)),'k')

figure(101)
cla
plot(t, imag(chirp_sig),'g')
hold on;

phi=@(t) 2*pi*(f0*t+(f1-f0)/(2*tau).*t.^2);
% 
% chirp = @(t,phi) ((t<(T-tau)/2 | t>(T+tau)/2).*0 ...
%                     +(t> (T-tau)/2 & t<=(T+tau)/2).*exp(1j.*phi));
% 
% phi=@(t,r) 2*pi*(f0*t+(f1-f0)/(2*tau).*t.^2);

phi_r=@(t,r) -(2*pi/c)*(f0+(f1+f0)/(2*tau).*t).*r;


chirp = @(t,phi) ((t<(T-tau)/2 | t>(T+tau)/2).*0 ...
                    +(t> (T-tau)/2 & t<=(T+tau)/2).*exp(1j.*phi));


chirp_r = @(t,phi,phi_r) ((t<(T-tau)/2 | t>(T+tau)/2).*0 ...
                    +(t> (T-tau)/2 & t<=(T+tau)/2).*exp(1j.*phi)).*exp(1j.*phi_r);
                
r=6;       
cla;
hold on;

sig=chirp(t,phi(t));

plot(t, imag(chirp_r(t-r/c,phi(t), phi_r(t-r/c,r))),'b')
plot(t, imag(chirp_r(t-r/c,phi(t), phi_r(t,r))),'g')
sig_reverse=conj(chirp(t-r/c,phi(t-r/c)));

plot(t, imag(chirp(t-r/c,phi(t))),'b')
plot(t, imag(conj(chirp(t-r/c,phi(t-r/c)))),'r')

plot(t, abs(chirp(t-r/c,phi(t-r/c)).*conj(chirp(t-r/c,phi(t-r/c)))),'k')

figure(102)
       
apollo=fft(real(chirp(t,phi(t))));
plot(0:1/T:44.1E3/2-1/T, apollo(1:end/2),'r')

