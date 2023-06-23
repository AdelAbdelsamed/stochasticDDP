%% compute_J function computes the cost for a given state and control input trajectory 
% 
% Inputs:
% x         : Cell array containing the state trajectory
% u         : Cell array containing the control input trajectory
% x_des     : Vector of the desired terminal state
% cost      : Instance of the class Cost 
% dt        : Scalar indicating the discretization step 
%
% Outputs:
% J         : Scaladr indicating Cost J for a given state and control input trajectory 
function J = compute_J(x ,u ,x_des ,cost , dt)
    
    J = cost.phi_N(x{end}, x_des); % Compute Terminal cost 
    for k = 1:length(x)-1
        J = J + cost.ell(x{k}, u{k}, dt); % Running Cost
    end

end