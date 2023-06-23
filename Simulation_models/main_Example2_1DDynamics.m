%% Main Scirpt for 1D example model 2
clc
clear
close all

% Define the initial state
x_0 = [0.0];

% Define the desired target state 
x_des = [2.0];

% Define the time horizon
t_f = 3.0; % in sec

dt = 0.02;
% Define the no. of discretization steps
N = floor(t_f ./ dt);

% Maximum magnitude of control
u_max = [1.0]; % in rad/sec

% Initialize dynamics
fprintf("initializing determinstic dynamics of example model...\n")
doi = 0.05; % degree of instability
sigma_squared = 0.0; % Noise variance
dyn = Example2_1DDynamics(doi, sigma_squared);

% Initialize cost
fprintf("initializing quadratic cost function...\n")
Q_f = 10.0;
R = 1e-3;
cost = QuadraticCost_Fcn(Q_f, R);

% Number of SDDP iterations
num_iter = 500;

% Step size control parameter
alpha = 1;
% De-/Activate line search algorithm
line_search_activated = true;
% De-/Activate Feedforward clamping 
ff_clamping_activated = false;


%% Run the SDDP algorithm
fprintf('running sddp...\n')

tic;
opt_sol = SDDP(x_0, x_des, t_f, N, dyn, cost, u_max, num_iter, alpha, line_search_activated, ff_clamping_activated);
time_elapsed = toc;

fprintf('The SDDP algorithm took %d seconds.\n', time_elapsed)


if opt_sol.error ~= 0
    fprintf("Error: %s \n",opt_sol.error_type);
end

% Extract solution
t_opt = opt_sol.t; % discrete time stamp vector

% Obtain optimal state and control input trajectories
x_opt = zeros(length(opt_sol.x),length(opt_sol.x{1}));
u_opt = zeros(length(opt_sol.u),length(opt_sol.u{1}));
for i = 1:length(x_opt)
    x_opt(i,:) = opt_sol.x{i};
    u_opt(i,:) = opt_sol.u{i};
end


%% Simulate system with zero noise
zero_noise_activated = true;
system_traj = simulate_system(x_0, opt_sol.u, t_f, N, dyn, zero_noise_activated);

% Plot evolution of objective functions over number of iterations
figure;
title('Cost History')
semilogy(1:num_iter,opt_sol.J,'r','LineWidth',1)
grid on
xlabel('No. of iterations [-]')
ylabel('Cost')

% Plot the optimal state and control input trajectories
figure;
subplot(2,1,1)
plot(system_traj.t_sim, system_traj.x_sim,'Color',[0,0.7,0.9], 'LineWidth',1)
grid on
ylabel('State x1 [-]')
xlabel('Time [sec]')

subplot(2,1,2)
plot(system_traj.t_sim, system_traj.u_sim,'b--','LineWidth',1)
grid on
ylabel('Optimal control input [-]')
xlabel('Time [sec]')
sgtitle('Optimal Trajectories');



