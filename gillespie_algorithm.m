% The Gillespie Algorithm (with two populations, each with 4 types of reactions)
% Kevin Roberts
% January 2025

clear all
close all
clc

k1v = 0.001; % (1/s) rate
k2v = 0.01; % (1/s) rate
k3v = 1.2; % (1/s) rate
k4v = 1; % (1/s) rate

number_of_tests = 5;
legend_entries_1 = cell(1, number_of_tests + 1); % Preallocate legend entries
legend_entries_2 = cell(1, number_of_tests + 1); % Preallocate legend entries

% the time vector is discretized with steps of dt = 1. Inside of the loop,
time_vec = 1:100;
A_sum_vec = zeros(1, 100);
B_sum_vec = zeros(1, 100);

figure(1)
hold on
title("Population A over time")
ylabel("Population");
xlabel("Time in ms");

figure(2)
hold on
title("Population B over time")
ylabel("Population");
xlabel("Time in ms");


for j = 1:number_of_tests
    
    % initial populations
    A_initial = 0; 
    A_vec(1) = A_initial;
    B_initial = 0; 
    B_vec(1) = B_initial;
    
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
        alpha_2 = A_vec(i)*B_vec(i)*k2v;
        alpha_3 = k3v;
        alpha_4 = k4v;
        alpha_0 = alpha_1 + alpha_2 + alpha_3 + alpha_4;

        % Compute the time when the next chemical reaction takes place as t + tau
        tau = 1/alpha_0 * log(1/r1);

        % adding to the sum vector of A and B
        if (current_time + tau) > discretized_time
            A_sum_vec(discretized_time) = A_sum_vec(discretized_time) + A_vec(i);
            B_sum_vec(discretized_time) = B_sum_vec(discretized_time) + B_vec(i);
            discretized_time = discretized_time + 1;
        end
        
        % Compute the number of molecules at time t + tau for A population
        if (0 <= r2) && (r2 < alpha_1/alpha_0)
            A_vec(i + 1) = A_vec(i) - 2;
        elseif (alpha_1/alpha_0 <= r2) && (r2 < (alpha_1 + alpha_2)/alpha_0)
            A_vec(i + 1) = A_vec(i) - 1;
        elseif ((alpha_1 + alpha_2)/alpha_0 <= r2) && (r2 < (alpha_1 + alpha_2 + alpha_3)/alpha_0)
            A_vec(i + 1) = A_vec(i) + 1;
        else % when ((alpha_1 + alpha_2 + alpha_3)/alpha_0 <= r2 < 1)
            A_vec(i + 1) = A_vec(i);
        end
        
        % Compute the number of molecules at time t + tau for A population
        if (0 <= r2) && (r2 < alpha_1/alpha_0)
            B_vec(i + 1) = B_vec(i);
        elseif (alpha_1/alpha_0 <= r2) && (r2 < (alpha_1 + alpha_2)/alpha_0)
            B_vec(i + 1) = B_vec(i) - 1;
        elseif ((alpha_1 + alpha_2)/alpha_0 <= r2) && (r2 < (alpha_1 + alpha_2 + alpha_3)/alpha_0)
            B_vec(i + 1) = B_vec(i);
        else % when ((alpha_1 + alpha_2 + alpha_3)/alpha_0 <= r2 < 1)
            B_vec(i + 1) = B_vec(i) + 1;
        end
        
        t_vec(i + 1) = t_vec(i) + tau;
        
        current_time = current_time + tau;
        i = i + 1;
    end
        
    % Plotting
    figure(1)
    stairs(t_vec, A_vec);
    hold on;
    legend_entries_1{j} = sprintf('Simulation %d', j); % Create legend entry
    
    figure(2)
    stairs(t_vec, B_vec);
    hold on;
    legend_entries_2{j} = sprintf('Simulation %d', j); % Create legend entry

end

avg_vec1 = A_sum_vec./number_of_tests;
avg_vec2 = B_sum_vec./number_of_tests;

figure(1)
plot(time_vec, avg_vec1, 'LineWidth', 2, 'Color', 'k')
legend_entries_1{end} = 'Mean of populations';
legend(legend_entries_1, "Interpreter", "latex");

figure(2)
plot(time_vec, avg_vec2, 'LineWidth', 2, 'Color', 'k')
legend_entries_2{end} = 'Mean of populations';
legend(legend_entries_2, "Interpreter", "latex");
