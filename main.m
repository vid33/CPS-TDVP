


clear;

timerval = zeros(1, 1);

h_list = 0.5;

for zzz=1:numel(h_list);

XXZ = false;  %standard pauli convention == xxz; flipped conventions (i.e. sz=[01;10] sx=[10;0-1] correspond to zzx)

IMAGINARY_TIME = true; %otherwise real-time

PRECONDITION = true; %Gram matrix inv w/ preconditioner or without?

TEXT_INFO = true; %dump info to command window
GRAPHIC_INFO = false; %live plot of energy/tangent vector norm..

bicg_maxiter =10000; bicg_tolerance = 1e-9;

siteOverlap = 2;
D = 2^siteOverlap;

%quantum ising

if XXZ == true;
    sx = [0,1;1,0];
    sz = [1,0;0,-1];
elseif XXZ == false %zzx
    sz = [0,1;1,0];
    sx = [1,0;0,-1];
end

J = 1;

h = h_list(zzz);
%h = 1;

if IMAGINARY_TIME == false %real-time ev. specify t_max at which to stop
   t_max = 10; 
elseif IMAGINARY_TIME == true %imag-time ev. value for norm^2 of tangent vector at which we stop imaginary t evolution
    tolerance = 1e-20;
    if XXZ == true
        fOut = sprintf('data/Ising_xxz/sites=%d/sites=%dJ=%dh=%d.mat',siteOverlap, siteOverlap, J, h); %save gs data to file
    elseif XXZ == false
        fOut = sprintf('data/Ising_zzx/sites=%d/sites=%dJ=%dh=%d.mat',siteOverlap, siteOverlap, J, h); %save gs data to file
    end
end
    

if exist('C', 'var') == 0;
    C = cpxrand(D,D); 
end

C_init = C;
t = 0; currentStep = 1;
dt = 0.2;



tic;
mainLoop_TDVP;
save(fOut);
toc;
%timerval(zz) = toc;  
end
