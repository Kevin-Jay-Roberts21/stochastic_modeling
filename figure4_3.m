% Reproducing figure 4.3
% Kevin Roberts
% April 2025

clear all
close all
clc

D = 0.0001; % (mm^2/sec) diffusion constant
h = 0.025; % (mm) space step

d = D/h^2; % (1/s) same as D/h^2
K = 40; % loop count (number of compartments)
X_0 = 0.4; % creating a random value to start with (could be anything)

% running the time loop
for i = 1:K
    
    if (i <= K/2)
        X = X_0 - h/2;
    else 
        X = X_0 + h/2;
    end
    
    total_time = 0;
    k = 1;
    A_plot(i, k) = X; % (in mm)  columns represent time, rows represent corresponding position for a certain A
    time_plot(i, k) = total_time; % (in sec) columns represent time, rows represent discrete time for a certain A
    
    while (k < 320) % 320 is arbitrary total time, could go for longer
        
        % creating random numbers
        r1 = rand;
        r2 = rand;
        
        % defining alpha_0 and tau
        alpha_0 = 2 * d;
        tau = 1/alpha_0 * log(1/r1);
        
        % updating the time
        total_time = total_time + tau;
        
        
        if (r2*alpha_0 < d) % move to the right
            X = X - h;
        else % move to the left
            X = X + h;
        end
        
        % applying boundary situations
        if (X < 0)
            X = h/2; 
        end
        
        if (X > 1)
            X = 1 - h/2;
        end
        
        % updating k and A_plot, time_plot vectors
        k = k + 1;
        A_plot(i, k) = X;
        time_plot(i, k) = total_time;
    end
end

% converting the time plot to be in minutes
time_plot=time_plot/60;

% plotting all the data
figure(1);
for i = 1:K
    plot = stairs(A_plot(i,:), time_plot(i,:));
    hold on
end

xlabel('$x$ [mm]','interpreter','latex');
ylabel('time [min]','interpreter','latex');
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'YTick',[0 2 4 6 8 10 12 14]);
axis([0 1 0 14]);
    
    
    
    
    
    
    