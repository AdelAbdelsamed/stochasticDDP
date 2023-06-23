%% This script defines the class Parafoil Dynamics which is a subclass of the abstract class Dynamics
classdef QuadraticCost_Fcn < Cost_Fcn


    properties
        Q_f; % Terminal State Penalty Weights
        R;   % Control Input Weights
    end


    methods
       
        % Constructor initializing objects of class QuadraticCost_Fcn
        function obj = QuadraticCost_Fcn(Q_f, R)
            obj.Q_f = Q_f;
            obj.R = R;
        end

        % Terminal Cost
        function phi_N = phi_N(obj, x_f, x_d)
            phi_N = 0.5 .*(x_f - x_d).' *obj.Q_f * (x_f - x_d);
        end

        % Derivatives of terminal state cost
        function phi_N_x = phi_N_x(obj, x_f, x_d)
            phi_N_x = obj.Q_f* (x_f - x_d);
        end

        function phi_N_xx = phi_N_xx(obj, ~, ~)
            phi_N_xx = obj.Q_f;
        end

        % Running Cost
        function ell = ell(obj, ~, u, dt)
            ell = (0.5 .* u' * obj.R * u) .* dt;
        end

        % Derivatives of the running cost
        function ell_x = ell_x(~, x, ~, ~)
            ell_x = zeros(length(x),1);
        end

        function ell_u = ell_u(obj, ~, u, dt)
            ell_u = (obj.R * u).* dt;
        end

        function ell_xx = ell_xx(~, x, ~, ~)
            ell_xx = zeros(length(x));
        end

        function ell_xu = ell_xu(~, x, u, ~)
            ell_xu = zeros(length(x),length(u));
        end

        function ell_ux = ell_ux(~, x, u, ~)
            ell_ux = zeros(length(u),length(x));
        end

         function ell_uu = ell_uu(obj, ~, ~, dt)
            ell_uu = obj.R .* dt;
         end


    end

end