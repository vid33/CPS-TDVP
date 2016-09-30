clear;

IMAGINARY_TIME = true; %otherwise real-time

TEXT_INFO = false; %print info in command window
GRAPHIC_INFO = false; %live plot of energy/tangent vector norm..

PRECONDITION = true;

bicg_maxiter =5000; bicg_tolerance = 1e-4;

siteOverlap = 4;
D = 2^siteOverlap;


phi = linspace(0,2*pi, 20);
theta = linspace(0, pi, 20);

%%%%%TESTING ANGLE CONVENTIONS
theta_x = pi/2; phi_x = 0;
sx = [cos(theta_x), sin(theta_x)*exp(-1i*phi_x);
    sin(theta_x)*exp(1i*phi_x), -cos(theta_x) ];

theta_y = pi/2; phi_y = pi/2;
sy = [cos(theta_y), sin(theta_y)*exp(-1i*phi_y);
    sin(theta_y)*exp(1i*phi_y), -cos(theta_y) ];

theta_z = 0; phi_z = rand;
sz = [cos(theta_z), sin(theta_z)*exp(-1i*phi_z);
    sin(theta_z)*exp(1i*phi_z), -cos(theta_z) ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi = linspace(0,2*pi, 20);
theta = linspace(0, pi, 20);

J = 1;
h = 1;

sx_zero = [0,1;1,0];
sz_zero = [1,0;0,-1];

%%% rotation about y axis
rot = linspace(-pi/2,pi/2, 200);
energy_scan  =zeros(1, numel(rot));

rot_power = 0; %2^(rot_power+1) copies of rot.

for kk =1: rot_power
   rot = [rot rot]; 
end


EE_info = zeros(1, numel(rot));


for kk=1:numel(rot)
    
    fprintf('At kk %d out of %d. \n', kk, numel(rot));
    
    Urot = [cos(rot(kk)/2), sin(rot(kk)/2);
            -sin(rot(kk)/2), cos(rot(kk)/2) ] ;
        
     sx = Urot*sx_zero*Urot';
     sz = Urot*sz_zero*Urot';


    if IMAGINARY_TIME == false %real-time ev. specify t_max at which to stop
       t_max = 10; 
    elseif IMAGINARY_TIME == true %imag-time ev. value for norm^2 of tangent vector at which we stop imaginary t evolution
        tolerance = 1e-20;
        fOut = sprintf('data/Ising_xxz/sites=%d/sites=%dJ=%dh=%d.mat',siteOverlap, siteOverlap, J, h); %save gs data to file
    end

    C = cpxrand(D,D); 
 
    t = 0; currentStep = 1;
    dt = 0.1;

    mainLoop_TDVP;

    energy_scan(kk) = energy;
    EE_info(kk) = min(eig(C.*conj(C)));

end

figure; plot(rot, energy_scan, 'x');

break;

%scan over theta

phi = pi/2;
for kk=1:numel(theta)

    sx = [cos(theta(kk)+(pi/2)), sin(theta(kk)+(pi/2))*exp(-1i*phi);
        sin(theta(kk)+(pi/2))*exp(1i*phi), -cos(theta(kk)+(pi/2)) ];

    sz = [cos(theta(kk)), sin(theta(kk))*exp(-1i*phi);
        sin(theta(kk))*exp(1i*phi), -cos(theta(kk)) ];


    if IMAGINARY_TIME == false %real-time ev. specify t_max at which to stop
       t_max = 10; 
    elseif IMAGINARY_TIME == true %imag-time ev. value for norm^2 of tangent vector at which we stop imaginary t evolution
        tolerance = 1e-20;
        fOut = sprintf('data/Ising_xxz/sites=%d/sites=%dJ=%dh=%d.mat',siteOverlap, siteOverlap, J, h); %save gs data to file
    end

    C = cpxrand(D,D); 

    t = 0; currentStep = 1;
    dt = 0.1;

    mainLoop_TDVP;

    energy_scan(kk) = energy;

end
  
