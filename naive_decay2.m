% Nonuniform time decay  (finding the time until the next reaction)
% Kevin Roberts 
% January 2025

close all
clear all
clc

k = 5; % (1/s) death rate
total_time = 1; % (s) total time
n = 101; % number of tau time steps (plus 1 because we don't count the first)
A_initial = 50; % (initial population) 

A_vec = zeros(1, n);
A_vec(1) = A_initial;
time_vec = zeros(1, n);

for i = 1:n
    r = rand;
    tau = 1/(A_vec(i)*k)*log(1/r);
    A_vec(i+1) = A_vec(i) - 1;
    time_vec(i+1) = time_vec(i) + tau;
end

figure(1)
plot(time_vec, A_vec);
ylim([0, 50])

ylabel("Population");
xlabel("Time in ms");
legend("Population $A(t)$", "Interpreter","latex");