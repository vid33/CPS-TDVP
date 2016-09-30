clear;

PLOT_MAX_TT_EIG = false; %check how for from zero max TT eigenvalue is when using VARY_Q

XXZ = false;  %standard pauli convention == xxz; flipped conventions (i.e. sz=[01;10] sx=[10;0-1] correspond to zzx)

t_min = 0;
t_max = 40;
delta_t = 0.01;
t_sol = t_min: delta_t: t_max;

siteOverlap = 2;
D = 2^siteOverlap;

MAX_PLOT_EIGV = D; %plot this many eigenvalues of TTmixed

%values of initial parameters
J_initial = 1;
h_initial = 0.75;

%quench to
J=1;
h=0.1;

quenchOverlap = zeros(1, numel(t_sol)-1);

%%%%%%%%%%%%%%%%%%%%%%%
if siteOverlap > 1
    twoPointCorr = zeros(1, numel(t_sol)-1);
    PP = 3; %sz  - - -(P-2) - - -sz correlator
    sz=[1 0;0 -1]; Sz = (1/2)*sz;
    sx=[0 1;1 0]; Sx = (1/2)*sx;
    twoPointLocal  = calculate_twoPoint_local( siteOverlap, Sx, PP);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eigEEmax = zeros(1, numel(t_sol)-1);

cpx_quenchOverlap = zeros(1, numel(t_sol)-1);

eigvalues = numel(t_sol-1, MAX_PLOT_EIGV); 

energy_quench = zeros(1, numel(t_sol));

for zz=1:numel(t_sol)

    fprintf('t_sol = %d\n', t_sol(zz));
    
    if XXZ == true
    fIn = sprintf('data_quench/Ising_xxz/sites=%d/sites=%dJ_initial=%dh_initial=%dJ=%dh=%dt=%d.mat', ...
                        siteOverlap, siteOverlap, J_initial, h_initial, J, h, t_sol(zz));
    
    elseif XXZ == false
    fIn = sprintf('data_quench/Ising_zzx/sites=%d/sites=%dJ_initial=%dh_initial=%dJ=%dh=%dt=%d.mat', ...
                        siteOverlap, siteOverlap, J_initial, h_initial, J, h, t_sol(zz));        
    end
    
    load(fIn, 'C', 'L', 'R', 'Ltrunc', 'Rtrunc', 'Bzero', 'energy');
    
    energy_quench(zz) = energy;
    
    if t_sol(zz) ==0
        C_init = C;
    end
    
    
     EEmixed = C_init.*conj(C);
          
     EE = C.*conj(C);
     
     %eigEEmax(zz) = max(real(eig(EE)));
     
    eigvalues_tmp = eig(EEmixed); [eigvalues_real_tmp, eigv_index]  = sort(eigvalues_tmp.*conj(eigvalues_tmp), 'descend');
    for xx=1:MAX_PLOT_EIGV
        eigvalues(zz,xx) = eigvalues_tmp(eigv_index(xx)); 
    end

    eigvalues_tmp = eig(EEmixed);
    [tmpsort, tmpindx] = sort(eigvalues_tmp.*conj(eigvalues_tmp), 'descend');
    quenchOverlap(zz) = eigvalues_tmp(tmpindx(1)).*conj(eigvalues_tmp(tmpindx(1)));
    cpx_quenchOverlap(zz) = eigvalues_tmp(tmpindx(1));

    if siteOverlap >1
        twoPointCorr(zz) = calculateTwoPointCorr_eff(C, L, R, twoPointLocal, PP);
    end
    %imag_quenchOverlap(zz) = imag(eigvalues_tmp(tmpindx(2)));
    
end


%figure; plott_sol, eigEEmax);
rate_function = -log(quenchOverlap)/siteOverlap;

figure; plot(t_sol/siteOverlap, rate_function);

%figure; plot(t_sol, log(cpx_quenchOverlap).*conj(log(cpx_quenchOverlap)));

figure; plot(t_sol/siteOverlap,(eigvalues.*conj(eigvalues)));

figure; plot(t_sol/siteOverlap, energy_quench);

if siteOverlap > 1
    figure; plot(t_sol/siteOverlap, twoPointCorr);
end
%figure; plot(t_sol, imag(eigvalues), 'x');



