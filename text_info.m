   

fprintf('********************************************\n');
fprintf('At step %d, at time %d\n', currentStep, t);
fprintf('Site overlap=%d, and D=%d\n', siteOverlap, D);
fprintf('energy is %d \n', energy);
fprintf('correlation length is %d \n', -1/log(Eeigvals(2)));
if IMAGINARY_TIME == true
fprintf('Norm of the tangent vector is %d (aim %d)\n', normDeltaC, tolerance); end
fprintf('delta energy is %d\n', deltaEnergy);
if t>0
fprintf('Ratio delta energy linear, delta energy %d\n', deltaEnergy/deltaEnergy_linear); end
fprintf('********************************************\n');