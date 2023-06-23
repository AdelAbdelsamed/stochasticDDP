%% Backward pass function
% Implementation of the backward pass 
% 
% Inputs:
% cost                 : Object of class cost_fcn representing the objective function of the problem
% dyn                  : Object of class dynamics representing the dynamics of the problem
% x                    : Cell array (Nx1) containing the states
% x_des                : Vector of the desired terminal state
% u                    : Cell array (Nx1) containing the control input
% dt                   : Scalar indicating the discretization interval 
% N                    : Scalar indicating the no. of discretization intervals
% kau                  : Scalar indicating approximation order of the dynamics
% V                    : Cell array (Nx1) containing the value function
% V_x                  : Cell array (Nx1) containing the first derivative of the value function
% V_xx                 : Cell array (Nx1) containing the second derivative of the value function
% mu                   : Scalar indicating the regularization parameter
% Q_x                  : Cell array (Nx1) containing the derivative of Q w.r.t to x
% Q_u                  : Cell array (Nx1) containing the derivative of Q w.r.t to u
% Q_xx                 : Cell array (Nx1) containing the derivative of Q w.r.t to xx
% Q_xu                 : Cell array (Nx1) containing the derivative of Q w.r.t to xu
% Q_ux                 : Cell array (Nx1) containing the derivative of Q w.r.t to ux
% Q_uu                 : Cell array (Nx1) containing the derivative of Q w.r.t to uu
%
% Outputs:
% V                    : Cell array (Nx1) containing the value function
% V_x                  : Cell array (Nx1) containing the first derivative of the value function
% V_xx                 : Cell array (Nx1) containing the second derivative of the value function
% k_ff                 : Cell array (Nx1) containing the feedforward gain
% K_fb                 : Cell array (Nx1) containing the feedback gain
% Q_x                  : Cell array (Nx1) containing the derivative of Q w.r.t to x
% Q_u                  : Cell array (Nx1) containing the derivative of Q w.r.t to u
% Q_xx                 : Cell array (Nx1) containing the derivative of Q w.r.t to xx
% Q_xu                 : Cell array (Nx1) containing the derivative of Q w.r.t to xu
% Q_ux                 : Cell array (Nx1) containing the derivative of Q w.r.t to ux
% Q_uu                 : Cell array (Nx1) containing the derivative of Q w.r.t to uu
% status_backward_pass : Scalar indicating the status of the backward pass
function [V, V_x, V_xx,k_ff ,K_fb, Q_x, Q_u, Q_xx, Q_xu, Q_ux, Q_uu, status_backward_pass] ...
                                          = backward_pass(cost, dyn, x, x_des, u, dt, N, kau, V, V_x, V_xx, mu, Q_x, Q_u, Q_xx, Q_xu, Q_ux, Q_uu)
    % Define temporary variables
    Q_xx_temp = cell(N,1);
    Q_ux_temp = cell(N,1);
    Q_xu_temp = cell(N,1);
    Q_uu_temp = cell(N,1);

    % Compute terminal value function and derivatives
    V(N) = cost.phi_N(x{N}, x_des);
    V_x{N} = cost.phi_N_x(x{N}, x_des);
    V_xx{N} = cost.phi_N_xx(x{N}, x_des);

    % Backward pass to obtain optimal control 
       for k = N-1:-1:1
            % Compute state-action value function derivatives
            Q_x{k} = cost.ell_x(x{k}, u{k}, dt) + dyn.A_t(x{k}, u{k}, dt).' * V_x{k+1} ...
                   + dyn.math_S_tilde(x{k}, u{k}, dt, V_xx{k+1});
            Q_u{k} = cost.ell_u(x{k}, u{k}, dt) + dyn.B_t(x{k}, u{k}, dt).' * V_x{k+1} ...
                   + dyn.math_U_tilde(x{k}, u{k}, dt, V_xx{k+1});
            Q_xx_temp{k} = cost.ell_xx(x{k}, u{k}, dt) + dyn.A_t(x{k}, u{k}, dt).' * V_xx{k+1}*dyn.A_t(x{k}, u{k}, dt) ...
                    + kau.*dyn.math_F(x{k}, u{k}, V_x{k+1}) + dyn.math_F_tilde(x{k}, u{k}, dt, V_xx{k+1}) ...
                    + dyn.math_M_tilde(x{k}, u{k}, dt, V_xx{k+1});
            Q_uu_temp{k} = cost.ell_uu(x{k}, u{k}, dt) + dyn.B_t(x{k}, u{k}, dt).' *V_xx{k+1}*dyn.B_t(x{k}, u{k}, dt) ...
                    + kau.*dyn.math_Z(x{k}, u{k}, V_x{k+1}) + dyn.math_Z_tilde(x{k}, u{k}, dt, V_xx{k+1}) ...
                    + dyn.math_G_tilde(x{k}, u{k}, dt, V_xx{k+1});
            Q_xu_temp{k} = cost.ell_xu(x{k}, u{k}, dt) + dyn.A_t(x{k}, u{k}, dt).' * V_xx{k+1}*dyn.B_t(x{k}, u{k}, dt) ...
                    + kau.*dyn.math_L(x{k}, u{k}, V_x{k+1}) + dyn.math_L_tilde(x{k}, u{k}, dt, V_xx{k+1}) ...
                    + dyn. math_N_tilde(x{k}, u{k}, dt, V_xx{k+1});
            Q_ux_temp{k} = Q_xu_temp{k}.';

            % Compute regularized matrices 
            Q_xx{k} = Q_xx_temp{k} + mu .* dyn.A_t(x{k}, u{k}, dt).' * eye(length(x{k})) * dyn.A_t(x{k}, u{k}, dt); % Add quadratic state cost
            Q_uu{k} = Q_uu_temp{k} + mu .* dyn.B_t(x{k}, u{k}, dt).' * eye(length(u{k})) * dyn.B_t(x{k}, u{k}, dt); % Add quadratic control cost
            Q_ux{k} = Q_ux_temp{k} + mu .* dyn.B_t(x{k}, u{k}, dt).' * eye(length(x{k})) * dyn.A_t(x{k}, u{k}, dt); % Add quadratic state/control cost
            Q_xu{k} = Q_ux{k}';


            % Check whether regularized Q_uu is positive definite
            if(~ all( eig(Q_uu{k}) > 0 ) )
               status_backward_pass = 1; % Indicate that backward pass has failed
               return
            end

           % Check whether Q_uu isnan
            if isnan(Q_uu{k})
               status_backward_pass = 2; % Indicate that backward pass has failed
               return
            end
        
           % Construct Gain Matrices:
           k_ff{k} = -inv(Q_uu{k})*Q_u{k};
           K_fb{k} = -inv(Q_uu{k})*Q_ux{k};

           % Compute the value function and its derivatives using the
           % modified value approximation
           V(k) = V(k) + 0.5.* k_ff{k}.' * Q_uu_temp{k} * k_ff{k} + k_ff{k}.' * Q_u{k}; 
           V_x{k} = Q_x{k} + K_fb{k}.' * Q_uu_temp{k} * k_ff{k} + K_fb{k}.' * Q_u{k} + Q_ux_temp{k}.' * k_ff{k};
           V_xx{k} = Q_xx{k} + K_fb{k}.' * Q_uu_temp{k} * K_fb{k} + K_fb{k}.' * Q_ux_temp{k} + Q_ux_temp{k}.' * K_fb{k} ;

       end
    
       status_backward_pass = 0;

end