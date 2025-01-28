% Multiple Reactions 2
% Kevin Roberts
% January 2025

clear all
close all
clc

k1v = 0.005; % (1/s) rate
k2v = 1; % (1/s) rate

number_of_tests = 5;
legend_entries = cell(1, number_of_tests + 1); % Preallocate legend entries

% the time vector is discretized with steps of dt = 1. Inside of the loop,
time_vec = 1:100;
sum_vec = zeros(1, 100);

figure(1)

for j = 1:number_of_tests
    
    A_initial = 0; % (initial population) 
    A_vec(1) = A_initial;
    t_vec(1) = 0;

    current_time = 0;
    i = 1;
    discretized_time = 1;
    
    while current_time < 100
        
        % Generate random variables
        r1 = rand;
        r2 = rand;

        % Compute propensity functions
        alpha_1 = A_vec(i)*(A_vec(i)-1) * k1v;
        alpha_2 = k2v;
        alpha_0 = alpha_1 + alpha_2;

        % Compute the time when the next chemical reaction takes place as t + tau
        tau = 1/alpha_0 * log(1/r1);

        % adding to the sum vector
        if (current_time + tau) > discretized_time
            sum_vec(discretized_time) = sum_vec(discretized_time) + A_vec(i);
            discretized_time = discretized_time + 1;
        end
        
        % Compute the number of molecules at time t + tau
        if r2 < alpha_1/alpha_0
            A_vec(i + 1) = A_vec(i) - 2;
        else
            A_vec(i + 1) = A_vec(i) + 1;
        end
        
        t_vec(i + 1) = t_vec(i) + tau;
        
        current_time = current_time + tau;
        i = i + 1;
    end
        
    % Plotting
    figure(1)
    stairs(t_vec, A_vec);
    hold on;
    legend_entries{j} = sprintf('Simulation %d', j); % Create legend entry

end

avg_vec = sum_vec./number_of_tests;

plot(time_vec, avg_vec, 'LineWidth', 2, 'Color', 'k')
legend_entries{end} = 'Mean of populations';


ylabel("Population");
xlabel("Time in ms");
legend(legend_entries, "Interpreter", "latex");


%%% Plotting the stationary distribution phi(n) function %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Getting the highest value from the average and rounding up (used to plot)

% max_value = ceil(max(sum_vec));
% 
% phi_vec = zeros(1, max_value); % initializing the vector
% 
% % populating the phi vector (n is the populations)
% for n = 1:max_value
%     phi_vec(n) = 1/factorial(n) * (k2v/k1v)^n * exp(-k2v/k1v);
% end
% 
% figure(2)
% pop_vec = 1:max_value;
% plot(pop_vec, phi_vec)
% ylabel("Stationary Distribution");
% xlabel("Population");
% legend('master equation');