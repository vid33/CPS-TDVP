    clear;

    J = 1;
    h = 1;

    XXZ = true;  %standard pauli convention == xxz; flipped conventions (i.e. sz=[01;10] sx=[10;0-1] correspond to zzx)

    if XXZ == false
        siteOverlap_list = 2:8;
    elseif XXZ == true
        siteOverlap_list = 2:10;
    end
    
    energy_list = zeros(1, numel(siteOverlap_list));
    corr_length = zeros(1, numel(siteOverlap_list));
    corr_length2 = zeros(1, numel(siteOverlap_list));
    entropy = zeros(1, numel(siteOverlap_list));
    ratio1 = zeros(1, numel(siteOverlap_list));
    
    Deff = zeros(1, numel(siteOverlap_list));
    
    for kk = 1:numel(siteOverlap_list)
   
        D = 2^siteOverlap_list(kk);
    
        if XXZ == true
            fIn = sprintf('data/Ising_xxz/sites=%d/sites=%dJ=%dh=%d.mat',siteOverlap_list(kk), siteOverlap_list(kk), J, h); %save gs data to file
        elseif XXZ == false
            fIn = sprintf('data/Ising_zzx/sites=%d/sites=%dJ=%dh=%d.mat',siteOverlap_list(kk), siteOverlap_list(kk), J, h); %save gs data to file
        end
    
        load(fIn, 'energy', 'C', 'L', 'R');
        energy_list(kk) = energy;
        
        EE = C.*conj(C);      
        eigenvalues = sort(real(eig(EE)));
        eigenvalues = flipud(eigenvalues);
    
        T_eig2 = eigenvalues(2); 
       T_eig3 = eigenvalues(3);
        corr_length(kk) = (-1/log(T_eig2))*siteOverlap_list(kk);
        corr_length2(kk) = (-1/log(T_eig3))*siteOverlap_list(kk);
        
        ratio1(kk) = (1/corr_length2(kk))/(1/corr_length(kk));
        
       %%%%%%%%%%%%%%%%%%%%%%%%
       rhor = C*diag(R(1,:))*C';
       rhol = diag(L(1,:));
       
       density_mat = rhol^(1/2)*rhor^(1/2);
       
       schmidt = svd(density_mat);
       schmidt_sq = schmidt.^2;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%      
%        schmidt_sq = L(1,:).*R(1,:);
        
        
        entropy(kk) = sum(-schmidt_sq.*log(schmidt_sq));
        
        Deff(kk) = (2^siteOverlap_list(kk));
        
    end
    
    exactEnergy_list = -1.2732395*ones(1, numel(siteOverlap_list));
    
    figure; plot(siteOverlap_list, energy_list, 'x-');

    %%%%%%%%%%%%%%%%
    figure; plot((siteOverlap_list), 1./(corr_length), 'x-');
    hold all;
    plot((siteOverlap_list), 1./(corr_length2), 'x-');
    hold off;
    %%%%%%%%%%%%%%%%%%%
    
    figure; plot(siteOverlap_list, entropy, 'x-');
        figure; plot(Deff, 1./corr_length, 'x-');
        
        figure; plot(log(corr_length), entropy, 'x-');
        
   deriv_entropy_vs_log_corr = zeros(1, numel(energy_list)-1); %derivative of entroly vs. log(mu); mu = correlation lenght
   
   for kk=1:numel(deriv_entropy_vs_log_corr)
       diffx = log(corr_length(kk+1)) - log(corr_length(kk));
       diffy = entropy(kk+1) - entropy(kk);
       
       deriv_entropy_vs_log_corr(kk) = diffy/diffx;
   end
        
    figure; plot(1./siteOverlap_list(2:end), deriv_entropy_vs_log_corr, 'x');
  
    
    figure; plot(siteOverlap_list, ratio1, 'x-');
    break;
        
%schmidt spectrum

spectrum = -log(sort(real(schmidt_sq)));
spectrum = flip(spectrum);

figure; plot(spectrum, 'x');


  %      figure; plot(Deff, ratio1, 'x-');

        
    
    
  % hold all;
  %  plot(siteOverlap_list, exactEnergy_list);
  %  hold off;
    
    
    
    
    