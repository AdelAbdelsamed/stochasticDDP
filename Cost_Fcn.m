%% This Script defines an abstract class for the cost function of the optimal control problem
classdef (Abstract) Cost_Fcn
    methods
        % Terminal Cost
        phi_N = phi_N(obj, x_f, x_d);
        % Derivatives of terminal cost
        phi_N_x = phi_N_x(obj, x_f, x_d);
        phi_N_xx = phi_N_xx(obj, x_f, x_d);

        % Running Cost
        ell = ell(obj,x, u, dt);
        % Derivatives of running cost
        ell_x = ell_x(obj,x, u, dt);
        ell_u = ell_u(obj,x, u, dt);
        ell_xx = ell_xx(obj,x, u, dt);
        ell_xu = ell_xu(obj,x, u, dt);
        ell_ux = ell_ux(obj,x, u, dt);
        ell_uu = ell_uu(obj,x, u, dt);

    end
end
