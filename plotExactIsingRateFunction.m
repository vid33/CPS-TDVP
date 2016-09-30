
%Evaluate and plot the exact Ising return functions using h_init and h_final as
%loaded


t_exact = linspace(0, 11, 1100);

rate_function_exact = zeros(1, numel(t_exact));

for kk=1:numel(t_exact)
    
    %rateFunctionIntegrand(k, t, J, h_init, h)
    rate_function_exact(kk) = integral(@(k)rateFunctionIntegrand(k, t_exact(kk), J, h_initial, h), 0, pi);   
    %rate_function(kk) = integral(@(k)rateFunctionIntegrand(k, t_exact(kk), J, h_initial, h), 0, pi, 'RelTol',0,'AbsTol',1e-12 );   
    
end

figure; plot(t_exact, 2*real(rate_function_exact));