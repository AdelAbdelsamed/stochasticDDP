%% Main Scirpt for Parafoil Homing Trajectory OptimizationS
clc
clear
close all
%% Set up the optimization problem
% Define the initial state
x_0 = [0.0 ; 0.0 ; 0.0; 0.0; 0.0; 0.0; 0.0];

% Define the desired target state 
x_des = [0.0; 100.0; 0.0; 0.0; 0.0; 0.0; 0.0];

% Define the time horizon
t_f = 5.0; % in sec

% Define the no. of discretization steps
dt = 0.02;
N = floor(t_f ./dt);

% Maximum magnitude of control
u_max = [5]; % in rad/sec

% Initialize cost
fprintf("initializing quadratic cost function...\n")
Q_f = [ 1.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0 0.0 0.0 0.0; ...
         0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; ...
         0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0 0.0];
R = 0.1;
cost = QuadraticCost_Fcn(Q_f, R);

% Number of SDDP iterations
num_iter = 5000;

% Step size control parameter
alpha = 1.0;
% De-/Activate line search algorithm
line_search_activated = true;
% De-/Activate Clamping of Control Input
ff_clamping_activated = false;
%% Generate Dryden Filter (Wind turbulence model)
% Filter parameters
sigma_turb = 0.1* 45; % Turbulence intensity
L = 400; % Turbulence length (Altitiude in feet)
sigma_wind = 0.1; % Covariance of white noise
x0_wind = sigma_wind.* randn(2,1); % Initial conditions
v = 15; % Horizontal velocity in m/s

% Create Laplace Variable
s = tf('s');

% Create Filter Transfer Function
DF_tf = sigma_turb * sqrt(L/ (pi * v)) * (1 + sqrt(3) * (L * s/ v) )/ ( (1 + (L * s/ v)^2) );
[DF_num, DF_den] = tfdata(DF_tf); % Obtain nom, den of filter

% Obtain ss matrices
[DF_A, DF_B, DF_C, DF_D] = tf2ss(DF_num{1},DF_den{1});

% Create minimal ss representation
DF_sys = minreal(ss(DF_A, DF_B, DF_C, DF_D));

% Generate multiple wind profiles
N_wind_profiles = 50; % No. of different realizations
y_wind = cell(N_wind_profiles);
% Plot wind profiles
figure; hold all
for i = 1:N_wind_profiles
    [y_wind{i}, t_wind] = generate_wind_profiles(x0_wind, N, dt, DF_A, DF_B, DF_C, DF_D, sigma_wind);
    plot(t_wind(1:end-1), y_wind{i}(1:end-1),'k', 'LineWidth',1);
end
% Highlight single wind profile
plot(t_wind(1:end-1), y_wind{N_wind_profiles}(1:end-1),'b', 'LineWidth',1.5)
grid on;
title('Wind profiles')
xlabel('Time [sec]')
ylabel('Wind intensity [m/sec]')

% Initialize dynamics for optimization problem
fprintf("initializing determinstic parafoil dynamics...\n")
v = 15; % Horizontal velocity in m/sec
sigma_squared = 1.0; % Variance
dyn = ParafoilDynamics(v, sigma_squared, DF_A, DF_B, DF_C );
%% Run the SDDP algorithm
fprintf('running sddp...\n')
tic;
opt_sol = SDDP4(x_0, x_des, t_f, N, dyn, cost, u_max, num_iter, alpha, line_search_activated, ff_clamping_activated);
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
zero_noise_activated = true;
system_traj = simulate_system(x_0, opt_sol.u, t_f, N, dyn, zero_noise_activated);

% Preproces sthe angle rate
theta_deg = rad2deg(system_traj.x_sim(3,:)); % Convert rad to deg
theta_deg_wrapped = mod(theta_deg,360); % Wrap angle around 360

% Obtain the windgusts 
wx = DF_C(1,1) .* system_traj.x_sim(4,:) + DF_C(1,2) .* system_traj.x_sim(5,:); % Wind gusts in x-direction
wy = DF_C(1,1) .* system_traj.x_sim(6,:) + DF_C(1,2) .* system_traj.x_sim(7,:); % Wind gusts in y-direction

% Plot evolution of objective functions over number of iterations
figure;
semilogy(1:num_iter,opt_sol.J,'r','LineWidth',1)
grid on
xlabel('No. of iterations [-]')
ylabel('Cost')
title('Cost History')

% Plot the optimal state and control input trajectories
figure;
subplot(2,3,1)
plot(system_traj.t_sim, system_traj.x_sim(1,:),'Color',[0,0.7,0.9], 'LineWidth',1)
hold on
plot(system_traj.t_sim, system_traj.x_sim(2,:),'Color',[0,0.2,0.9], 'LineWidth',1)
grid on
ylabel('Position in xy-plane [m]')
xlabel('Time [sec]')
legend('State x1', 'State x2')

subplot(2,3,2)
plot(system_traj.t_sim, wx ,'Color',[0.8500 0.3250 0.0980], 'LineWidth',1)
hold on
plot(system_traj.t_sim, wy,'Color',[0.4660 0.6740 0.1880], 'LineWidth',1)
ylabel('Wind intenisty [m/s]')
xlabel('Time [sec]')
legend('wind intensity in x-direction','wind intensity in y-direction')

subplot(2,3,3)
plot(system_traj.t_sim, system_traj.x_sim(3,:) ,'Color',[0,0.7,0.9], 'LineWidth',1)
grid on
ylabel('Theta [deg]')
xlabel('Time [sec]')

subplot(2,1,2)
plot(system_traj.t_sim, system_traj.u_sim,'b--','LineWidth',1)
grid on
ylabel('Optimal control input: Turn rate [rad/sec]')
xlabel('Time [sec]')
sgtitle('Optimal Trajectories');

