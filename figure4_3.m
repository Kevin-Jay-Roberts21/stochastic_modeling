% Reproducing figure 4.3
% Kevin Roberts
% April 2025

clear all
close all
clc

% creating the initial conditions and rates
kv1 = ;
kv2 = ;
kv3 = ;
kv4 = ;
d = 0.16; % (1/s)
K = 40; % number of compartments


A_1 = 0;
A_2 = 0;
A_3 = 0;
A_4 = 0;
A_5 = 0;
A_6 = 0;
A_7 = 0;
A_8 = 0;
A_9 = 0;
A_10 = 0;
A_11 = 0;
A_12 = 0;
A_13 = 0;
A_14 = 0;
A_15 = 0;
A_16 = 500;
A_17 = 500;

list_of_As = {A_1 , A_2, A_3, A_4, A_5, A_6, A_7, A_8, A_8, A_9, A_10, A_11, A_12, A_13, A_14, A_15, A_16, A_17}

A_s = zeros(1, 40);

A_s(16) = 500;
A_s(17) = 500;

% running the time loop
for k = 1:100
    r1 = rand;
    r2 = rand;
    
    alpha_0 = 2*d*K;
    
    tau = 1/alpha_0 * log(1/r1);
    
    sum = sum(A_s)/alpha_0;

    if r2 < sum
    
        
        
    else
        
    
    end
    
    
    
    
    % compute propensity function
    alpha_0 = 2 * d; 
    


end
