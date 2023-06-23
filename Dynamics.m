%% This Script defines an abstract class for the dynamics model of the system
classdef (Abstract) Dynamics

    methods
        % This function implements the deterministic dynamics
        f_deter = f(obj,x,u);
        % This function implements the stochastic dynamics
        F_stoch = F(obj,x ,u);

        
        % Discrete-Time Linearized Second-order model
        A_t = A_t(obj,x ,u ,dt);
        B_t = B_t(obj,x ,u ,dt);

        % Functions corresponding to the expectation of the exxpanded value function:
        % Terms arrising from the first-order term of the value function
        math_F_out = math_F(obj , x, u, V_x);
        math_Z_out = math_Z(obj , x, u, V_x);
        math_L_out = math_L(obj , x, u, V_x);

        % Terms arrising from the second-order term of the value function
        math_F_tilde_out = math_F_tilde(obj, x, u, dt, V_xx);
        math_L_tilde_out = math_L_tilde(obj, x, u, dt, V_xx);
        math_Z_tilde_out = math_Z_tilde(obj, x, u, dt, V_xx);
        math_U_tilde_out = math_U_tilde(obj, x, u, dt, V_xx);
        math_S_tilde_out = math_S_tilde(obj, x, u, dt, V_xx);

        math_M_tilde_out = math_M_tilde(obj, x, u, dt, V_xx);
        math_N_tilde_out = math_N_tilde(obj, x, u, dt, V_xx);
        math_G_tilde_out = math_G_tilde(obj, x, u, dt, V_xx);
        % Function returns the dimension of the noise
        getDimensionm = getDimensionm(obj);

    end
end