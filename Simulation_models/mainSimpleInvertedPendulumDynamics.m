%% Main Scirpt for Parafoil Homing Trajectory Optimization
clc
clear
close all
%%
% Define the initial state
x_0 = [-pi ; 0.0];

% Define the desired target state 
x_des = [0.0; 0.0];

% Define the time horizon
t_f = 4.0; % in sec

% Define the no. of discretization steps
dt = 0.01;
N = floor(t_f ./dt);

% Maximum magnitude of control
u_max = [5.0]; % in rad/sec

% Initialize dynamics
fprintf("initializing inverted pendulum dynamics...\n")
beta = 0.04; % Measurement parameter noise
sigma_squared = 1.0; % Noise Variance
dyn = SimpleInvertedPendulumDynamics(beta, sigma_squared);

% Initialize cost
fprintf("initializing quadratic cost function...\n")
Q_f = [10.0 0.0;0.0 0.0];
R = 1e-2;
cost = QuadraticCost_Fcn(Q_f, R);

% Number of SDDP iterations
num_iter = 10000;

% Step size control parameter
alpha = 1.0;
% De-/Activate line search algorithm
line_search_activated = true;
% De-/Activate Clamping of Control Input
ff_clamping_activated = false;


%% Run the SDDP algorithm
fprintf('running sddp...\n')

tic;
opt_sol = SDDP5(x_0, x_des, t_f, N, dyn, cost, u_max, num_iter, alpha, line_search_activated, ff_clamping_activated);
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
for i = 1:N
    x_opt(i,:) = opt_sol.x{i};
    u_opt(i,:) = opt_sol.u{i};
end


%% Simulate system with zero noise
% Plot evolution of objective functions over number of iterations
grayColor = [.7 .7 .7];
figure;
semilogy(1:num_iter,opt_sol.J,'r','LineWidth',1)
grid on
xlabel('No. of iterations [-]')
ylabel('Cost')
title('Cost History')

% Plot multiple instances and the mean of the optimal state trajectory
zero_noise_activated = false;
figure;
for j = 1:20
    system_traj = simulate_system(x_0, opt_sol.u, t_f, N, dyn, zero_noise_activated);
    
    % Preproces sthe angle rate
    theta_deg = rad2deg(system_traj.x_sim(1,:)); % Convert rad to deg
    theta_deg_wrapped = mod(theta_deg,360); % Wrap angle around 360
    
    % Plot the optimal state and control input trajectories
    subplot(3,1,1)
    plot(system_traj.t_sim, system_traj.x_sim(1,:),'Color',grayColor, 'LineWidth',1)
    hold on
    grid on
    ylabel('Angle [rad]')
    xlabel('Time [sec]')

    subplot(3,1,2)
    plot(system_traj.t_sim, system_traj.x_sim(2,:),'Color',grayColor, 'LineWidth',1)
    hold on
    grid on
    ylabel('Turn Rate [rad/sec]')
    xlabel('Time [sec]')
    
    subplot(3,1,3)
    plot(system_traj.t_sim, system_traj.u_sim,'b','LineWidth',1)
    hold on
    grid on
    ylabel('Optimal control input: Turn rate [rad/sec]')
    xlabel('Time [sec]')
    sgtitle('Optimal Trajectories');
end

zero_noise_activated = true;
system_traj = simulate_system(x_0, opt_sol.u, t_f, N, dyn, zero_noise_activated);
    
% Preproces sthe angle rate
theta_deg = rad2deg(system_traj.x_sim(1,:)); % Convert rad to deg
theta_deg_wrapped = mod(theta_deg,360); % Wrap angle around 360

% Plot the optimal state and control input trajectories
subplot(3,1,1)
plot(system_traj.t_sim, system_traj.x_sim(1,:),'b', 'LineWidth',1)
hold on
grid on
ylabel('Angle [rad]')
xlabel('Time [sec]')

subplot(3,1,2)
plot(system_traj.t_sim, system_traj.x_sim(2,:),'b', 'LineWidth',1)
hold on
grid on
ylabel('Turn Rate [rad/sec]')
xlabel('Time [sec]')

subplot(3,1,3)
plot(system_traj.t_sim, system_traj.u_sim,'b','LineWidth',1)
hold on
grid on
ylabel('Optimal control input: Turn rate [rad/sec]')
xlabel('Time [sec]')
sgtitle('Optimal Trajectories');




