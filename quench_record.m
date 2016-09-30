%Run quench

%%%NB INCREASE TOLERANCE - LOWERED NOW FOR quenchOverlap = 8

clear;

global    Bzero_GLOBAL;
%Bzero_GLOBAL = cpxrand(D*(D-1), 1);

XXZ = false;  %standard orientation is xxz; flipped conventions (i.e. sz=[01;10] sx=[10;0-1] correspond to zzx)

STORED = true; %record to file - only option for now

IMAGINARY_TIME = false; %vs real-time evolution

siteOverlap = 8;
D = 2^siteOverlap;

%quantum ising
if XXZ == true;
    sx = [0,1;1,0];
    sz = [1,0;0,-1];
elseif XXZ == false %zzx
    sz = [0,1;1,0];
    sx = [1,0;0,-1];
end

bicg_maxiter=20000;
bicg_tolerance=5e-9;

options_quench = odeset('RelTol', 1e-10, 'AbsTol', 1e-13, 'Stats', 'on');

t_min =25.98;
t_max = 26;
delta_t = 0.01;
tspan_ode = [0 delta_t];

t_sol = t_min: delta_t: t_max;

%values of initial parameters
J_initial = 1;
h_initial = 1.5;

%quench to
J=1;
h=0.2;

if t_min == 0    
    fprintf('LOADING INITIAL\n');
    if XXZ == true 
        fIn = sprintf('data/Ising_xxz/sites=%d/sites=%dJ=%dh=%d.mat',siteOverlap, siteOverlap, J_initial, h_initial); 
    elseif XXZ == false
        fIn = sprintf('data/Ising_zzx/sites=%d/sites=%dJ=%dh=%d.mat',siteOverlap, siteOverlap, J_initial, h_initial);
    end
    load(fIn, 'C', 'L', 'R', 'Ltrunc', 'Rtrunc', 'Bzero');

else
    if XXZ == true
        fIn = sprintf('data_quench/Ising_xxz/sites=%d/sites=%dJ_initial=%dh_initial=%dJ=%dh=%dt=%d.mat', ...
                        siteOverlap, siteOverlap, J_initial, h_initial, J, h, t_min);
    elseif XXZ == false
        fIn = sprintf('data_quench/Ising_zzx/sites=%d/sites=%dJ_initial=%dh_initial=%dJ=%dh=%dt=%d.mat', ...
                siteOverlap, siteOverlap, J_initial, h_initial, J, h, t_min);
    end
    load(fIn, 'C', 'L', 'R', 'Ltrunc', 'Rtrunc', 'Bzero', 'Bzero_GLOBAL');
        
end

if exist('Bzero', 'var') == 0 
    Bzero_GLOBAL = cpxrand(D*(D-1), 1);
else
    Bzero_GLOBAL = Bzero;
end


[ Hlocal  ] = calculate_H_local( siteOverlap, sx, sz, J, h );
[ energy ] = calculateEnergy_eff(C, L, R, Hlocal);

if STORED == true   && t_min == 0
    if XXZ == true
        fOut = sprintf('data_quench/Ising_xxz/sites=%d/sites=%dJ_initial=%dh_initial=%dJ=%dh=%dt=%d.mat', ...
                    siteOverlap, siteOverlap, J_initial, h_initial, J, h, t_min);
    elseif XXZ == false
        fOut = sprintf('data_quench/Ising_zzx/sites=%d/sites=%dJ_initial=%dh_initial=%dJ=%dh=%dt=%d.mat', ...
                    siteOverlap, siteOverlap, J_initial, h_initial, J, h, t_min);        
    end
    save(fOut, 'C', 'L', 'R', 'Ltrunc', 'Rtrunc', 'Bzero', 'energy');
end


for kk=1:(numel(t_sol)-1)
    
    fprintf('Time is %d\n',  t_sol(kk+1));

    [ Hlocal  ] = calculate_H_local( siteOverlap, sx, sz, J, h );

    solution = ode113(@(t,x) calculateDeltaC_ODE(x, D, Hlocal, bicg_maxiter, bicg_tolerance), tspan_ode, ... 
        reshape(transpose(C), D^2, 1), options_quench );    
        
    Csol = solution.y;
    
    C = transpose(reshape(Csol(1:D^2, end), D, D));
      
    [ C, Eeigvals, L, R ] = calculateEigenvectors_inv(C);
    
    Rtrunc = R(2:end, :);
    Ltrunc = L(2:end, :);
    
    energy_old = energy;
   [ energy ] = calculateEnergy_eff(C, L, R, Hlocal);

   %diag(L*transpose(R)) 
   
   
    if STORED == true
        if XXZ == true
            fOut = sprintf('data_quench/Ising_xxz/sites=%d/sites=%dJ_initial=%dh_initial=%dJ=%dh=%dt=%d.mat', ...
                        siteOverlap, siteOverlap, J_initial, h_initial, J, h, t_sol(kk+1));
        elseif XXZ ==false
            fOut = sprintf('data_quench/Ising_zzx/sites=%d/sites=%dJ_initial=%dh_initial=%dJ=%dh=%dt=%d.mat', ...
                        siteOverlap, siteOverlap, J_initial, h_initial, J, h, t_sol(kk+1));
        end
       
        Bzero = Bzero_GLOBAL;
        save(fOut, 'C', 'L', 'R', 'Ltrunc', 'Rtrunc', 'Bzero_GLOBAL', 'energy', 'bicg_maxiter', 'bicg_tolerance');
    end
   
     deltaEnergy = energy - energy_old;       
     fprintf('deltaEnergy is %d\n', deltaEnergy);
 
end

