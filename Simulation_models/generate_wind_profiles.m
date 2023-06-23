%% Function generating wind profiles
%
% Input Arguments:
% x_0         : Initial position
% N           : Discretization points
% dt          : Time intervall
% A           : System Matrix
% B           : Input Matrix
% C           : Output Matrix
% D           : "Durchgriff"-Matrix
% sigma_wind  : Covariance of Gaussian white noise which serves as input
%
% Output Arguments:
% y_wind      : Cell array containg the wind profiles

function [y_wind, t_wind] = generate_wind_profiles(x_0, N, dt, A, B, C, D, sigma_wind)
    
    % Allocate matrices

    % Time vector
    t_wind = 0:dt:dt*(N-1);
    % State tarjectory 
    x = cell(N,1);
    x{1} = x_0;
    
    % Control input trajectory
    u = cell(N,1);
    u{N} = zeros(1,1);

    % Output trajectory 
    y_wind = zeros(N,1);

    % Generate the wind profiles
    for k = 1:N-1
        % Simulate system one step forward
        x{k+1} = (eye(numel(x_0)) + A.*dt) * x{k} + sqrt(dt).* B .* sigma_wind.* randn();
        % Obtain wind profile
        y_wind(k) = C * x{k};
    end

end