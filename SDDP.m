%% Stochastic Differential Dynamic Programming Algorithm by Theodouru
% Implementation of the SDDP Algorithm with backtracking line search and
% the first regularization method
% 
% Inputs:
% x_0                  : Vector representing the initial condition
% x_des                : Vector of the desired terminal state
% t_f                  : Scalar indicating the final time 
% N                    : Scalar indicating the no. of discretization intervals
% dyn                  : Object of class dynamics representing the dynamics of the problem
% cost                 : Object of class cost_fcn representing the objective function of the problem
% u_max                : Vector indicating the maximum admissibile control input for each input
% num_iter             : Scalar indicating the number of iterations of the SDDP algorithm
% alpha                : Scalar value for the step size control
% line_search_activated: Scalar value  true = Activates the line search algorithm, false for deactivation
% ff_clamping_activated: Scalar value  true = Activates the feed forward clamping, false for deactivation
% varargin             : Cell with one entriy:
%           varargin{1}: Initial control input trajectory u_init: defined as a matrix with dimension (N,p) with p inputs 
%
% Outputs:
% opt_sol   : structure containing the locally optimal state and control trajectory
function opt_sol = SDDP(x_0,x_des,t_f, N, dyn, cost, u_max,num_iter, alpha, line_search_activated, ff_clamping_activated, varargin)    
    %% Allocate arrays for algorithm
    % Obtain time array
    t = linspace(0.0,t_f, N);
    dt = t(2) - t(1); % Discretization interval

    % Cost history 
    J = zeros(num_iter,1);
    J_temp = 0; % Cost of trajectory generated during forward pass before line search condition is checked 

    % State Trajectory
    x = cell(N, 1);
    x_plus = cell(N, 1); % New trajectory
    x{1} = x_0; % Initial Condition
    x_plus{1} = x_0;

    % Control Input Trajectory
    u = cell(N, 1);

    % Value function and its derivatives
    V = zeros(N,1);
    V_x = cell(N,1);
    V_xx = cell(N,1);

    % State-action value function derivatives
    Q_x = cell(N, 1);
    Q_u = cell(N, 1);
    Q_xx = cell(N, 1);
    Q_uu = cell(N, 1);
    Q_xu = cell(N, 1);
    Q_ux = cell(N, 1);
    
    % Helper variables used:
    % For temporary control inputs:
    u_temp = cell(N,1); % Control input trajectory generated before being accepted by the line search algorithm
    u_temp{N} = zeros(size(u_max));
    % For temporary state trajectories
    x_plus_temp = cell(N,1); % State trajectory generated before being accepted by the line search algorithm
    x_plus_temp{1} = x_0;
    % For temporary state-action value functions to ensure the positive definiteness of the right lower block [Qxx Qxu;Qux Quu] 
    Q_uu_temp = cell(N,1); 
    Q_xx_temp = cell(N,1);

    
    %% Tweak the SDDP parameters
    % SDDP Parameters
    roh = 0.7;             % Backtracking decay
    B = 40;                % Maximum no. of backtracking iterations 
    kau = 1;               % Determines the order of the expanded dynamics: 0 for 1st order, 1 for 2nd order
    tolFun = 1e-9;         % Termination criterion: Min. reduction in objective function
    tolGradNorm = 1e-7;    % Termination criterion: Min. of the norm of the Gradient 

    %% Initialize the DDP with a random input sequence

    fprintf("initializing random control sequence...\n")
    
    % Check to see if an initial control input is provided, if not then randomly generate one using u_max 
    if (nargin < 12)
        % Generate random control sequence
        for k = 1:N-1
            u{k} = 2 .* u_max .* rand(length(u_max),1) - u_max;
        end
        u{N} = zeros(numel(u_max), 1); % Enforce u_N = 0 
    else
        % Assign given control input
        u_init = varargin{1};
        for k = 1:N
            u{k} = u_init(k,:)';
        end
    end
    
    fprintf("generating initial trajectory...\n")


    % Use the random control sequence to obtain initial trajectory
    for k = 1:N-1
          x_plus{k+1} = x_plus{k} + dyn.f(x_plus{k},u{k}).*dt; % Instead of computing the mean of multiple realizations we just zero out the noise
    end
    
    fprintf("initial trajectory obtained ...\n")

    %% Perform the SDDP algorithm

    fprintf("launching sddp algorithm...\n")

    for i = 1:num_iter
        fprintf("SDDP iteration %d out of %d ...\n",i,num_iter)
        
        if ( i > 1)
            backtracking_iter = 1; % No. of backtracking iterations
            alpha_dec = alpha; % Initialize alpha with one

            % Repeat forward pass until line search condition is satisfied
             while (true)
                % Forward Pass
                for k = 1:N-1
                    % Compute control update feed-forward and feed-back
                    du_ff = -inv(Q_uu{k}) * Q_u{k};
                    du_fb = -inv(Q_uu{k}) * Q_ux{k} * (x_plus_temp{k} - x{k});
    
                    % Limit feed forward control with naive clamping from
                    % [Control-Limited Differential Dynamic Programming, Tassa et.al]
                    if (ff_clamping_activated) 
                        du_ff = min(u_max, max(-u_max,du_ff + u{k})) - u{k};
                    end
    
                    % Update control (Step-Size Control)
                    u_temp{k} = u{k} + alpha_dec .* (du_ff) + du_fb; 
                    
                    % Compute next state in trajectory with new control  (Generate new trajectory)  
                    x_plus_temp{k+1} = x_plus_temp{k} + dyn.f(x_plus_temp{k},u_temp{k}).*dt;
    
                    % Return error if new trajectory is not defined
                    if isnan(x_plus_temp{k+1})
                      opt_sol = generate_solution_struct(x, u, t, J, Q_u, Q_uu, Q_ux, 1);
                      return
                    end

                end

                % Backtracking line search 
                J_temp = compute_J(x_plus_temp ,u_temp ,x_des ,cost , dt); % Compute the cost for the newly generated trajectory
                        
                % Check if cost decreases or maximum no. of back iterations is reached or line search is activated 
                if( J_temp < J(i-1) || backtracking_iter > B || ~line_search_activated )
                    % Accept updated trajectories
                     u = u_temp;
                     x_plus = x_plus_temp;
                     break;
                else % If cost is not decreasing or maximum no. of backtracking iterations not achieved, then reduce the step size parameter by the backtracking decay factor
                    alpha_dec = alpha_dec .* roh;
                    backtracking_iter = backtracking_iter + 1; 
                end

             end

        end

        % Update current trajectory
        x = x_plus;

        % Compute total cost for iteration i
        J(i) = compute_J(x ,u ,x_des ,cost , dt);
        fprintf('The cost at stage %d is: %d \n',i,J(i));
        
         % Check whether SDDP converged
        if (i > 1)
            if ( abs(J(i) - J(i-1)) < tolFun )
                fprintf('SDDP algorithm converged: dV < dV_min, preparing results...\n')
                opt_sol = generate_solution_struct(x, u, t, J, Q_u, Q_uu, Q_ux, 0);
                return;
            end
        end

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
                Q_xu{k} = cost.ell_xu(x{k}, u{k}, dt) + dyn.A_t(x{k}, u{k}, dt).' * V_xx{k+1}*dyn.B_t(x{k}, u{k}, dt) ...
                        + kau.*dyn.math_L(x{k}, u{k}, V_x{k+1}) + dyn.math_L_tilde(x{k}, u{k}, dt, V_xx{k+1}) ...
                        + dyn. math_N_tilde(x{k}, u{k}, dt, V_xx{k+1});
                Q_ux{k} = Q_xu{k}.';
                
                % Check whether Q_uu is NaN
                if isnan(Q_uu_temp{k})
                  opt_sol = generate_solution_struct(x, u, t, J, Q_u, Q_uu, Q_ux, 2);
                  return
                end

                % Initialize Q_uu with mu = 0
               Q_uu{k} = Q_uu_temp{k};
               Q_xx{k} = Q_xx_temp{k};
               mu = 1.0; % Regularization parameter

               % Ensure the matrix H = [Qxx Qxu;Qux Quu] is positive definite by adding quadratic 
               % state and control costs to Q_xx and Q_uu stagewise
               H = [Q_xx{k} Q_xu{k};Q_ux{k} Q_uu{k}]; % Construct matrix H
               % Check the matrix [Quu] for positive definiteness
               while( ~all(eig(H) > 0) )
                   Q_uu{k} = Q_uu_temp{k} + mu .* eye(size(u{k})); % Add control cost
                   Q_xx{k} = Q_xx_temp{k} + mu .* eye(size(x{k})); % Add control cost
                   mu = 1.2 .* mu; % Increase mu if, not positive definite

                   % Check whether mu is infinity and regularization becomes infeasible (diverging costs)
                   if ( isinf(mu) )
                       opt_sol = generate_solution_struct(x, u, t, J, Q_u, Q_uu, Q_ux, 3);
                       return;
                   end

                   H = [Q_xx{k} Q_xu{k};Q_ux{k} Q_uu{k}]; % Reconstruct modified H matrix
               end

                % Compute the value function and its derivatives
                V(k) = V(k) - 0.5.* Q_u{k} * (Q_uu{k} \ Q_u{k}); 
                V_x{k} = Q_x{k} - Q_xu{k} * (Q_uu{k} \ Q_u{k});
                V_xx{k} = Q_xx{k} - Q_xu{k} * (Q_uu{k} \ Q_ux{k});
           end

           % Compute gradient norm Q_u
           Qu_norm = 0;
           for z = 1:N-1
              Qu_norm = Qu_norm + norm(Q_u{z});
           end
           
           % Terminate ?: Gradient norm is sufficiently small
           if( abs(Qu_norm) < tolGradNorm)
                fprintf('SDDP algorithm converged: Qu_norm < Qu_norm_min, preparing results...\n')
                opt_sol = generate_solution_struct(x, u, t, J, Q_u, Q_uu, Q_ux, 0);
                return;
           end

    
    end

    % SDDP Algorithm converged.. Generate solution structure
    fprintf('SDDP algorithm is done, preparing results...\n')
    
    opt_sol = generate_solution_struct(x, u, t, J, Q_u, Q_uu, Q_ux, 0);
end