%% Helper Function for SDDP Algorithm
% Generates solution structure
%
% Inputs
% x               : locally optimal state trajectory
% u               : locally optimal control sequence
% t               : discretized time stamps
% J               : iteration history of cost function
% Q_u, Q_uu, Q_ux : derivatives of state action value function
% error           : zero if no error, nonzero else
%
% Outputs
%
% sol : solution structure
function sol = generate_solution_struct(x, u, t, J, Q_u, Q_uu, Q_ux, error)

    % Solution structure
    sol = struct;
    
    sol.error = error;
    % Create error message
    if sol.error == 0
        sol.error_type = 'No Errors!';
    elseif sol.error == 1
        sol.error_type = 'Generated trajectory contains a NaN instance!';
    elseif sol.error == 2
        sol.error_type = 'Generated Q_uu is NaN!';
    elseif sol.error == 3 
        sol.error_type = 'Regularization parameter mu = inf! Diverging costs expected!';
    else
        sol.error_type = 'Error is not known!';
    end
    
    sol.x = x; 
    sol.u = u;
    sol.t = t;
    sol.dt = t(2) - t(1);
    sol.J = J;
    sol.Q_u = Q_u;
    sol.Q_uu = Q_uu;
    sol.Q_ux = Q_ux;
    
    return
end