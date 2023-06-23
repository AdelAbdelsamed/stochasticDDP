%% This script defines the class Inverted Pendulum Dynamics which is a subclass of the abstract class Dynamics
classdef SimpleInvertedPendulumDynamics < Dynamics
    properties (Access = public)
        beta; % Constant velocity in xy-plane
        sigma_squared; % Noise variance
    end

    properties (Constant)
       m = 1; % Dimension of noise 
    end

    methods
        
    % Constructor initializing objects of class SimpleInvertedPendulumDynamics
    function obj = SimpleInvertedPendulumDynamics(beta, sigma_squared)
        obj.beta = beta;
        obj.sigma_squared = sigma_squared;
    end

    % Discretized original general system

    % Function capturing the deterministic dynamics
    function f_deter = f(~,x,u)
        dx1 = x(2);
        dx2 = 4*sin(x(1)) + u;
        f_deter =  [dx1; dx2];
    end

    % Function capturing the stochastic dynamics
    function F_stoch = F(obj,~ ,u)
        F_stoch = [0 ; obj.beta* u];
    end

    % Discretized linearized second-order system
    function A_t = A_t(~,x ,~ ,dt)
        % Compute the gradient of the deterministic dynamics w.r.t x
        grad_f_x = zeros(length(x));
        grad_f_x(1,2) = 1;
        grad_f_x(2,1) = 4 * cos(x(1));

        A_t = eye(length(x)) + grad_f_x .* dt;
    end

    function B_t = B_t(~,x ,u ,dt)
        % Compute the gradient of the deterministic dynamics w.r.t u
        grad_f_u = zeros(length(x),length(u));
        grad_f_u(2,1) = 1;

        B_t = grad_f_u .* dt;
    end

    % Functions corresponding to the expectation of the expanded value function:
    % Terms arrising from the first-order term of the value function
    function math_F_out = math_F(~ , x, ~, V_x)
        math_F_out = zeros(length(x));
        % Specify the gradients
        grad_f_xx_i = cell(length(x),1); % Cell array containing the gradients grad_xx_f_i
        grad_f_xx_i{1} = zeros(length(x));
        grad_f_xx_i{2} = zeros(length(x));
        grad_f_xx_i{2}(1,1) = -4 * sin(x(1));
        
        for i = 1:length(x)
            math_F_out = math_F_out + grad_f_xx_i{i} .* V_x(i);
        end
    end

     function math_Z_out = math_Z(~ , x, u, V_x)
        math_Z_out = zeros(length(u));
        % Specify the gradients
        grad_f_uu_i = cell(length(x),1); % Cell array containing the gradient grad_uu_f_i
        grad_f_uu_i{1} = zeros(length(u));
        grad_f_uu_i{2} = zeros(length(u));
        
        for i = 1:length(x)
            math_Z_out = math_Z_out + grad_f_uu_i{i} .* V_x(i);
        end
     end

     function math_L_out = math_L(~ , x, u, V_x)
        math_L_out = zeros(length(x),length(u));
        % Specify the gradients
        grad_f_xu_i = cell(length(x),1); % Cell array containing the gradient grad_uu_f_i
        grad_f_xu_i{1} = zeros(length(x),length(u));
        grad_f_xu_i{2} = zeros(length(x),length(u));
        
        for i = 1:length(x)
            math_L_out = math_L_out + grad_f_xu_i{i}.*V_x(i);
        end
     end

     % Terms arrising from the second-order term of the value function
     function math_F_tilde_out = math_F_tilde(~,x, ~, ~, ~)
        math_F_tilde_out = zeros(length(x));
     end

     function math_L_tilde_out = math_L_tilde(~,x, u, ~, ~)
        math_L_tilde_out = zeros(length(x),length(u));
     end

     function math_Z_tilde_out = math_Z_tilde(obj,~, ~, dt, V_xx)
        math_Z_tilde_out = obj.sigma_squared .* dt .* [0 obj.beta] * V_xx * [0; obj.beta];
     end

     function math_U_tilde_out = math_U_tilde(obj, ~, u, dt, V_xx)
        math_U_tilde_out = obj.sigma_squared .* dt .* [0 obj.beta] * V_xx * [0; obj.beta * u];
     end

     function math_S_tilde_out = math_S_tilde(~,x, ~, ~, ~)
        math_S_tilde_out = zeros(length(x),1);
     end

     function math_M_tilde_out = math_M_tilde(~, x, ~, ~, ~)
        math_M_tilde_out = zeros(length(x));
     end

     function math_N_tilde_out = math_N_tilde(~, x, u, ~, ~)
        math_N_tilde_out = zeros(length(x),length(u));
     end

     function math_G_tilde_out = math_G_tilde(~, ~, u, ~, ~)
        math_G_tilde_out = zeros(length(u));
     end
    
     % Function returns the dimension of the noise
     function getDimensionm = getDimensionm(~)
        getDimensionm = SimpleInvertedPendulumDynamics.m;
     end



end

end