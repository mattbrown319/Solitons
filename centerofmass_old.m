%center of mass calculation


index = 1;
CenterOfMass = [];
numSteps = MAXTIME/del_t;
while index < numSteps
    OneCenterOfMass = trapz(x,x'.*(abs(uArray(:,index)).^2))/trapz(x,(abs(uArray(:,index))).^2); 
    CenterOfMass = [CenterOfMass OneCenterOfMass];
    index = index + 1;
end

plot(CenterOfMass)

max(abs(fft(CenterOfMass')))