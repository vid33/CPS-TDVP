deltaEnergy = 0;
normDeltaC = 1;

%Create a stop button to exit loop nicely
figh = figure;

global IS_ABORTED;
IS_ABORTED = false;

btn = uicontrol('style', 'pushb', 'string', 'Abort', ...
                'callback', @doAbort);
drawnow;

if GRAPHIC_INFO == true
    figureName = sprintf('Ising, J=%d, h=%d', J, h);
    energyPlot = figure('Name', figureName); 
    figure(energyPlot);
    set(0,'CurrentFigure', energyPlot);
    xlabel('time');
    ylabel('energy');
end

while 1

    [ C, Eeigvals, L, R ] = calculateEigenvectors_inv(C);
    
    Rtrunc = R(2:end, :);
    Ltrunc = L(2:end, :);
    
%    [ G ] = calculateG( C, R, L );
   % [ Gtest ] = calculateG2( C, R, L, Rtrunc, Ltrunc );
 %   [ Gecon ] = calculateGecon( C, R, L, Rtrunc, Ltrunc );
        
    [ Hlocal  ] = calculate_H_local( siteOverlap, sx, sz, J, h );
    [ Wvec ] = calculateDeltaC( C, L, R, Ltrunc, Rtrunc, Eeigvals, Hlocal);
    
    if currentStep == 1
        Bzero = cpxrand(D*(D-1), 1);
    end
    
    if PRECONDITION == true
        B = G_inv_action_precond( Wvec, C, L, R, Ltrunc, Rtrunc, bicg_maxiter, bicg_tolerance, Bzero );
    else
        B = G_inv_action( Wvec, C, L, R, Rtrunc, bicg_maxiter, bicg_tolerance, Bzero );
    end
    Bzero = B;

    
    Bmat = transpose(reshape(B, D, D-1));
    deltaC = (transpose(Rtrunc)*Bmat*L).*(1./conj(C));

    normDeltaC = B'*G_times_vec( B, C, L, R, Rtrunc );
       
    gradOverlap = B'*Wvec;
    deltaEnergy_linear =  -dt*(gradOverlap + conj(gradOverlap));
    
    if t >0
        energy_old = energy;
    end
    
    %energy = calculateEnergy(C, L(1,:), transpose(R(1,:)), J, h);
    [ energy ] = calculateEnergy_eff(C, L, R, Hlocal);
    
    if t>0
        deltaEnergy = energy - energy_old; 
    end
    
    if rem(currentStep,1)==0      
        if TEXT_INFO == true
           text_info; 
        end
    end
    if GRAPHIC_INFO == true
        if rem(currentStep, 250) == 0 || currentStep == 1;
            figure(energyPlot);
            %if rem(currentStep,100) == 0;
            %    clf(energyPlot);
            %    xlabel('time');
            %    ylabel('energy'); 
            %end
            hold on;
            plot(t, energy, 'x');
            drawnow;
            hold off;
        end
    end
        
    drawnow;
    if IS_ABORTED == true
        close(figh);
        break;
    end
    
    if IMAGINARY_TIME == true && deltaEnergy > 1e-12 && t > 0
        fprintf('Energy increasing.. stopping..\n');
        %close(figh);
       % break;
    end
    
    if IMAGINARY_TIME == true && normDeltaC < tolerance
        fprintf('C converged to desired tolerance. Saving data. \n');
        %save(fOut); 
        close(figh);
        break;
    end
    
    %backup in case things go wrong
    if rem(currentStep,100)==0
        save('data/tmp_save.mat');
    end 
    
    if IMAGINARY_TIME == true
        C = C - dt*deltaC;
    elseif IMAGINARY_TIME ==false
        C = C - 1i*dt*deltaC;
    end
    
    t = t + dt; currentStep = currentStep+1;
    
end

