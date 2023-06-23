%% This script defines the class Simple 1d example model which is a subclass of the abstract class Dynamics
%% This script is intended to test the implementation of the SDDP algorithm 
classdef Example2_1DDynamics < Dynamics

    properties (Constant)
       m = 1; % Dimension of noise 
    end

    properties
        alpha; % Degree of instability
        sigma_squared; % Noise variance
    end

    methods

    % Constructor
    function obj = Example2_1DDynamics(alpha, sigma_squared)
        obj.alpha = alpha;
        obj.sigma_squared = sigma_squared;
    end

    % Discretized original general system

    % Function capturing the deterministic dynamics
    function f_dyn = f(obj,x,u)
        dx = - obj.alpha .* x.^2 + u;
        f_dyn =  [dx];
    end

    % Function capturing the stochastic dynamics
    function F_stoc = F(~,x ,~)
        F_stoc = x.^2;
    end

    % Discretized linearized second-order system
    function A_t = A_t(obj,x ,~ ,dt)
        % Compute the gradient of the deterministic dynamics w.r.t x
        A_t = eye(length(x)) - obj.alpha .* 2 .* x .* dt;
    end

    function B_t = B_t(~,~ ,~ ,dt)
        % Compute the gradient of the deterministic dynamics w.r.t u
        B_t = dt;
    end

    % Functions corresponding to the expectation of the expanded value function:
    % Terms arrising from the first-order term of the value function
    function math_F_out = math_F(obj , ~, ~, V_x)
        math_F_out = obj.alpha .* 2 *V_x;

    end

     function math_Z_out = math_Z(~ , ~, u, ~)
        math_Z_out = zeros(length(u));   
     end

     function math_L_out = math_L(~ , x, u, ~)
        math_L_out = zeros(length(x),length(u));
     end

     % Terms arrising from the second-order term of the value function
     function math_F_tilde_out = math_F_tilde(obj ,x, ~, dt, V_xx)
        math_F_tilde_out = 4 .* obj.sigma_squared .* dt .* V_xx .* x.^2;
     end

     function math_L_tilde_out = math_L_tilde(~,x, u, ~, ~)
        math_L_tilde_out = zeros(length(x),length(u));
     end

     function math_Z_tilde_out = math_Z_tilde(~,~, u, ~, ~)
        math_Z_tilde_out = zeros(length(u));
     end

     function math_U_tilde_out = math_U_tilde(~,~, u, ~, ~)
        math_U_tilde_out = zeros(length(u),1);
     end

     function math_S_tilde_out = math_S_tilde(obj ,x ,~ , dt, V_xx)
        math_S_tilde_out = 2 .* obj.sigma_squared .* dt .* V_xx .* x.^3;
     end

     function math_M_tilde_out = math_M_tilde(obj ,x ,~ , dt, V_xx)
        math_M_tilde_out = 2 .* obj.sigma_squared .* dt .* V_xx .* x.^2;
     end

     function math_N_tilde_out = math_N_tilde(~, x, u, ~, ~)
        math_N_tilde_out = zeros(length(x),length(u));
     end

     function math_G_tilde_out = math_G_tilde(~, ~, u, ~, ~)
        math_G_tilde_out = zeros(length(u));
     end
    
     % Function returns the dimension of the noise
     function getDimensionm = getDimensionm(~)
        getDimensionm = DeterministicParafoilDynamics.m;
     end



end

end