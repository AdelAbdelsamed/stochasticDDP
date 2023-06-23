%% This script defines the class Parafoil Dynamics which is a subclass of the abstract class Dynamics
classdef ParafoilDynamics < Dynamics
    properties (Access = public)
        v;              % Constant velocity in xy-plane
        sigma_squared;  % Noise variance
        A_DrydenFilter; % System Matrix for Dryden Filter
        B_DrydenFilter; % Input Matrix for Dryden Filter
        C_DrydenFilter; % Output Matrix for Dryden Filter
    end

    properties (Constant)
       m = 2; % Dimension of noise 
    end

    methods
        
    % Constructor initializing objects of class DeterministicParafoilDynamics
    function obj = ParafoilDynamics(v, sigma_squared, A_DrydenFilter, B_DrydenFilter, C_DrydenFilter)
        obj.v = v;
        obj.sigma_squared = sigma_squared;
        obj.A_DrydenFilter = A_DrydenFilter;
        obj.B_DrydenFilter = B_DrydenFilter;
        obj.C_DrydenFilter = C_DrydenFilter;
    end

    % Discretized original general system

    % Function capturing the deterministic dynamics
    function f_deter = f(obj,x,u)
        dx1 = obj.v .* cos(x(3)) + obj.C_DrydenFilter*[x(4); x(5)];
        dx2 = obj.v .* sin(x(3)) + obj.C_DrydenFilter*[x(6); x(7)];
        dx3 = u;
        dx4 = obj.A_DrydenFilter(1,:)* [x(4); x(5)];
        dx5 = obj.A_DrydenFilter(2,:)* [x(4); x(5)];
        dx6 = obj.A_DrydenFilter(1,:)* [x(6); x(7)];
        dx7 = obj.A_DrydenFilter(2,:)* [x(6); x(7)];
        f_deter =  [dx1; dx2; dx3; dx4; dx5; dx6; dx7];
    end

    % Function capturing the stochastic dynamics
    function F_stoch = F(obj,~ ,~)
        F_stoch = [0 0; 0 0; 0 0; obj.B_DrydenFilter(1,:) 0; obj.B_DrydenFilter(2,:) 0;...
                     0 obj.B_DrydenFilter(1,:); 0 obj.B_DrydenFilter(2,:)];
    end

    % Discretized linearized second-order system
    function A_t = A_t(obj,x ,~ ,dt)
        % Compute the gradient of the deterministic dynamics w.r.t x
        grad_f_x = zeros(length(x));
        grad_f_x(1,3) = - obj.v .* sin(x(3));
        grad_f_x(1,4:5) = obj.C_DrydenFilter;
        grad_f_x(2,3) = obj.v .* cos(x(3));
        grad_f_x(2,6:7) = obj.C_DrydenFilter;
        grad_f_x(4:5,4:5) = obj.A_DrydenFilter;
        grad_f_x(6:7,6:7) = obj.A_DrydenFilter;

        A_t = eye(length(x)) + grad_f_x .* dt;
    end

    function B_t = B_t(~,x ,u ,dt)
        % Compute the gradient of the deterministic dynamics w.r.t u
        grad_f_u = zeros(length(x),length(u));
        grad_f_u(3,1) = 1;

        B_t = grad_f_u .* dt;
    end

    % Functions corresponding to the expectation of the expanded value function:
    % Terms arrising from the first-order term of the value function
    function math_F_out = math_F(obj , x, ~, V_x)
        math_F_out = zeros(length(x));

        math_F_out(1,3) = - obj.v .* cos(x(3)) .* V_x(1);
        math_F_out(2,3) = - obj.v .* sin(x(3)) .* V_x(2);
    end

     function math_Z_out = math_Z(~ , ~, u, ~)
        math_Z_out = zeros(length(u));
     end

     function math_L_out = math_L(~ , x, u, ~)
        math_L_out = zeros(length(x),length(u));
     end

     % Terms arrising from the second-order term of the value function
     function math_F_tilde_out = math_F_tilde(~,x, ~, ~, ~)
        math_F_tilde_out = zeros(length(x));
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
        getDimensionm = ParafoilDynamics.m;
     end



end

end