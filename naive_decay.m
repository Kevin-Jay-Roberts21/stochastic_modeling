% Stochastic Algorithm for Decay
% Kevin Roberts
% January 2025

clear all
close all
clc

k = 5; % (1/s) death rate
dt = 0.001; % (s) time step
total_time = 2; % (s) total time
n = total_time/dt + 1; % number of time steps (plus 1 because we don't count the first)
A_initial = 50; % (initial population) 

A_vec = zeros(1, n);
A_vec(1) = A_initial;

for i = 1:n-1
    r = rand;
    if r < A_vec(i)*k*dt
        A_vec(i+1) = A_vec(i) - 1;
    else
        A_vec(i+1) = A_vec(i);
    end
end

figure(1)
t = linspace(0, total_time, n);
plot(t, A_vec);
ylim([0, 50])

ylabel("Population");
xlabel("Time in ms");
legend("Population $A(t)$", "Interpreter","latex");
