%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%Extract static RF signal transmitted through acoustically modulated barrier
%
%Application of Dan Mark's algorithmic approach 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function Aprime = removeModulatedRFSignal(measurements,ModulationFunct,t,yy,zz)

dt=t(2)-t(1);
T=t(end)-t(1);





%calculate cross correlation

Gamma_yz0=zeros(length(zz),length(yy));
for nz = 1:length(zz)
    for ny= 1:length(yy)
%         Gamma_yz0 = xcorr(measurements(nz,ny,:),f_modulation(t));
        Gamma_yz0(nz,ny) = sum((squeeze(squeeze(measurements(nz,ny,:))).*ModulationFunct(t)'*dt));
    end
end

Ayz = mean(measurements,3);

alphayz  = Ayz.*Gamma_yz0;
betaxy   = Gamma_yz0.*conj(Gamma_yz0);

%return signal not attributed to modulated scatterers
Aprime=Ayz-alphayz./alphayz.*Gamma_yz0;