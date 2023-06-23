%% Function simulating a nonlinear stochastic system
% Generating a system trajectory using a control input sequence
% 
% Inputs:
% x_0                  : Vector representing the initial condition
% u                    : Cell array containing the control input sequence
% t_f                  : Scalar indicating the final time 
% N                    : Scalar indicating the no. of discretization intervals
% dyn                  : Object of class dynamics representing the dynamics of the problem
% zero_noise_activated : Boolean controlling the noise 
%
% Outputs:
% system_traj   : structure containing the locally optimal state and control trajectory as well as time vector
function system_traj = simulate_system(x_0, u, t_f, N, dyn, zero_noise_activated)
    
    % Simulated state trajectory 
    x_sim = zeros(length(x_0),N);
    x_sim(:,1) = x_0;

    % Simulated noise trajectory 
    xi_sim = zeros(dyn.getDimensionm,N);

    % Simulated control input trajectory 
    u_sim = zeros(length(u{1}),N);
    u_sim(:,N) = zeros(length(u{1}),1);

    % Create time vector 
    t_sim = linspace(0,t_f,N);
    dt = t_sim(2) - t_sim(1);

    % Noise variance
    if (zero_noise_activated)
        noise_variance = 0.0;
    else
        noise_variance = dyn.sigma_squared;
    end

    for i = 1:N-1
        % Obtain control input at i
        u_sim(:,i) = u{i};
        % Generate random variable xi_sim at i
        xi_sim(:,i) = sqrt(noise_variance) .* randn(dyn.getDimensionm,1);
        % Simulate one step forward i+1
        x_sim(:,i+1) = x_sim(:,i) + dyn.f(x_sim(:,i),u{i}).* dt + dyn.F(x_sim(:,i),u{i})* xi_sim(:,i) .*sqrt(dt);
    end
    
    % Create ouput struct
    system_traj.x_sim = x_sim;
    system_traj.xi_sim = xi_sim;
    system_traj.u_sim = u_sim;
    system_traj.t_sim = t_sim;

end