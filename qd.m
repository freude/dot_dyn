classdef qd
    % write a description of the class here.
    
    properties
        % define the properties of the class here, (like fields of a struct)
        num_vn;
        num_dots;
        C;
        Cv;
        Nmax;
        T;
        RS;
        RD;
    end
    
    methods
        % methods, including the constructor are defined in this block
        
        function obj = qd(C,Cdv,Nmax,T)
            % class constructor
            if (nargin == 1)
                
                obj.T=0.3;
                obj.num_vn=6;
                obj.num_dots=2;
                obj.Nmax=5;
                
                obj.Cv=-[Cg   Cg    100   0    100   Cg1;...
                          Cc   Cc     0   100   100   300].*q;                %cooupling between dots and voltage nodes

                obj.C=-[0.0   Ci;...
                         Ci   0.0].*q;                              %cooupling between dots
                     
            elseif (nargin==2)
                obj.T=0.3;
                s=size(Cdv);
                obj.num_dots=s(1);
                obj.num_vn=s(2);
                obj.Nmax=5;
                obj.Cv=Cdv;                %cooupling between dots and voltage nodes
                obj.C=C;                              %cooupling between dots                
            else
                s=size(Cdv);
                obj.num_dots=s(1);
                obj.num_vn=s(2);
                obj.T=T;
                obj.Nmax=Nmax;
                obj.Cv=Cdv;                %cooupling between dots and voltage nodes
                obj.C=C;
            end
        end
              
        function [E,varargout]=SD_SET(obj,V,N)
            
            q=1.6e-19;       
            
            mu_minus=zeros(1,obj.num_dots);
            mu_plus=zeros(1,obj.num_dots);
            Vd=zeros(1,obj.num_dots);
            
            E=0.5*(obj.C\(-q.*N-obj.Cv*V))'*(-q.*N-obj.Cv*V); %free energy
            
            for j=1:obj.num_dots
                Naux=N(j)-1;
                mu_minus(j)=E-0.5*(obj.C\(-q.*Naux-obj.Cv*V))'*(-q.*Naux-obj.Cv*V);
                Naux=N(j)+1;
                mu_plus(j)=0.5*(obj.C\(-q.*Naux-obj.Cv*V))'*(-q.*Naux-obj.Cv*V)-E;            
                Vd(j)=obj.C\(-q.*N-obj.Cv*V);
            end;
            
            varargout{1} = Vd;
            varargout{2} = mu_minus;
            varargout{3} = mu_plus;
            
        end
        
        function G=rates(obj,i,j,flag,Vsd)
            
            q=1.6e-19;
            
            if flag=='S'
                R=obj.RS;
            else
                R=obj.RD;
            end;
            
            if ((i-j)==1)
                [~,~,~,mu_plus]=obj.SD_SET(V,N);
                E=mu_plus-Vsd;
            elseif ((i-j)==-1)
                [~,~,mu_minus,~]=obj.SD_SET(V,N);
                E=Vsd-mu_minus;
            else
                E=q;
            end;
            
            G=1/q/R/q*ff(E);        
            
        end
        
        
    end
end


