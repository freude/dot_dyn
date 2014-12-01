classdef semicond_semicond
    % the class implements the quantum dot system
    
    properties
        % define the properties of the class here, (like fields of a struct)
        phi;
    end
    
    methods
        % methods, including the constructor are defined in this block
        
        function obj = semicond_semicond()
            q=1.6e-19;
            obj.phi=0.02*q; 
        end
        
        function T=transparency(obj,d,E,V)
            h=6.63e-34;
            m=0.19*9.1e-31;
            q=1.6e-19;
            ph=obj.phi*(1-0.25*q*V/obj.phi)^2;
            gamma=sqrt(2)*pi*pi/h*sqrt(m/obj.phi)*d;

            T=1/(1+exp(gamma*(ph-E)));
            %T=0.00001;
        end
    end
end