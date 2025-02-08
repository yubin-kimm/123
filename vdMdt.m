% Gallagher

function result = vdMdt(T,B)
C0=13.67388+2.718523*T+0.1278728*T^2-2.5072373*10^-3*T^3;
C1=-2.980799+0.9144735*T-7.3461875*10^-2*T^2+1.1088194*10^-3*T^3+9.8684995*10^-6*T^4;
C2=3.166695-0.6364962*T+4.8334930*10^-2*T^2-1.0840200*10^-3*T^3+4.3787923*10^-6*T^4;
C3=-0.3461423+9.0944469*10^-2*T-6.9979527*10^-3*T^2+1.3868879*10^-4*T^3;
C4=1.3560383*10^-2-4.1968897*10^-3*T+3.3840450*10^-4*T^2-6.8436516*10^-6*T^3;

dC0dT=2.718523+2*0.1278728*T-3*2.5072373*10^-3*T^2;
dC1dT=0.9144735-2*7.3461875*10^-2*T+3*1.1088194*10^-3*T^2+4*9.8684995*10^-6*T^3;
dC2dT=-0.6364962+2*4.8334930*10^-2*T-3*1.0840200*10^-3*T^2+4*4.3787923*10^-6*T^3;
dC3dT=9.0944469*10^-2-2*6.9979527*10^-3*T+3*1.3868879*10^-4*T^2;
dC4dT=-4.1968897*10^-3+2*3.3840450*10^-4*T-3*6.8436516*10^-6*T^2;

DEN=C0+C1*B+C2*B^2+C3*B^3+C4*B^4;
dDENdT=dC0dT+dC1dT*B+dC2dT*B^2+dC3dT*B^3+dC4dT*B^4;

% cs are [g] unit! So multiply 1000

result = -B*dDENdT/DEN^2*1000;
end







