classdef semicond_metal
    % the class implements the quantum dot system
    
    properties
        % define the properties of the class here, (like fields of a struct)
        phi;
        phi_s; % work_functions
        phi_m; % work_functions
        R;
    end
    
    methods
        % methods, including the constructor are defined in this block
        
        function obj = semicond_metal()
            q=1.6e-19;
            obj.phi_s=0*4.85*q+5.1*q; % work_functions
            obj.phi_m=5.1*q; % work_functions
            obj.phi=0.5*(obj.phi_m+obj.phi_s); % work_functions
            obj.R=0.1e-9;
        end
        
        function T=transparency(obj,d,E,V)
            h=1.05e-34;
            m=9.1e-31;
            q=1.6e-19;
            T=exp(-2*(d+obj.R)*2/3*sqrt(2*m/h/h)*((obj.phi_m-E+V*q)^(3/2)-(obj.phi_s-E)^(3/2))/(obj.phi_m-obj.phi_s+V*q));
            %T=0.00001;
        end
    end
end
