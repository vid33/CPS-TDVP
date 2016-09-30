%standard orientation is xxz; flipped conventions (i.e. sz = [01;10] correspond to zzx)

clear; 

XXZ = false;

siteOverlap = 2;

D = 2^siteOverlap;


if XXZ == true;
    sx = [0,1;1,0];
    sz = [1,0;0,-1];
elseif XXZ == false %zzx
    sz = [0,1;1,0];
    sx = [1,0;0,-1];
end
J = 1;
     
hscan =  0.9:0.01:1.3;
%hscan = 1.5;
sx_expectation = zeros(1, numel(hscan));
 
[ sxlocal  ] = calculate_sx_local( siteOverlap, sx );


 for kk=1:numel(hscan)
     
    if XXZ == true 
        fIn = sprintf('data/Ising_xxz/sites=%d/sites=%dJ=%dh=%d.mat',siteOverlap, siteOverlap, J, hscan(kk)); %save gs data to file
    elseif XXZ == false 
            fIn = sprintf('data/Ising_zzx/sites=%d/sites=%dJ=%dh=%d.mat',siteOverlap, siteOverlap, J, hscan(kk)); %save gs data to file
    end
    
    load(fIn, 'C', 'L', 'R');
 
    [ sx_expectation(kk) ] = calculateEnergy_eff(C, L, R, sxlocal);

 
 end

 
 figure; plot(hscan, abs(sx_expectation), 'x-');