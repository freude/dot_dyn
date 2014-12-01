classdef sqd
    % the class implements the quantum dot system
    
    properties
        % define the properties of the class here, (like fields of a struct)
        num_vn;
        num_dots;
        num_states;
        C;
        Cv;
        Nmax;
        T;
        RS;
        RD;
    end
    
    methods
        % methods, including the constructor are defined in this block
        
        function obj = sqd(C,Cdv,Nmax,T)
            obj.RS=2000000;
            obj.RD=2000000;
            obj.num_states=2;
            % class constructor
            if (nargin == 0)
                q=1.6e-19;
                obj.T=0.03;
                obj.num_vn=3;
                obj.num_dots=1;
                obj.Nmax=15;
                obj.Cv=-[300  200   200].*q;                %cooupling between dots and voltage nodes
                obj.C=300*q;                              %cooupling between dots                
            elseif (nargin==2)
                obj.T=0.3;
                s=size(Cdv);
                obj.num_dots=s(1);
                obj.num_vn=s(2);
                obj.Nmax=15;
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
            % function computes electrostatic properties of the system for
            % given gate voltages V and charge configuration N
            
            q=1.6e-19;       
            
            mu_minus=zeros(1,obj.num_dots);
            mu_plus=zeros(1,obj.num_dots);
            Vd=zeros(1,obj.num_dots);
            
            E=0.5*(obj.C\(-q.*N-obj.Cv*V))'*(-q.*N-obj.Cv*V); %free energy
            
            for j=1:obj.num_dots
                Naux=N(j)-1;
                mu_minus(j)=E-0.5*(obj.C\(-q.*Naux-obj.Cv*V))'*(-q.*Naux-obj.Cv*V);
                Naux=N(j)+1;
                mu_plus(j)=   0.5*(obj.C\(-q.*Naux-obj.Cv*V))'*(-q.*Naux-obj.Cv*V)-E;            
                Vd(j)=obj.C\(-q.*N-obj.Cv*V);
            end;
            

            varargout{1} = mu_minus;
            varargout{2} = mu_plus;
            varargout{3} = Vd;
            
        end
        
        function G=rates(obj,i,j,flag,V,i1,j1)
            % the function computes the transition rate from charge state i
            % to the charge state j at given voltages V
            % flag specifies what lead is involved in the transition
            
            q=1.6e-19;            
                        
            if flag=='S'
                R=obj.RS;
                Vsd=V(2);
            else
                R=obj.RD;
                Vsd=V(3);
            end;
            
            [~,mu_minus,mu_plus]=obj.SD_SET(V, i);
            
            if ((i-j)==1) %release an electron
                E=q*Vsd+mu_minus;
            elseif ((i-j)==-1) %take an electron
                E=q*Vsd-mu_plus;
            else
                E=q;
            end;
            
            if nargin>5 % there are several well defined energy states
                G(j1,j2)=1/q/R*ff_q(E,obj.T);
            else        % there is several well defined energy states
                G=1/q/R/q*ff(E,obj.T);
            end;
            
        end
        
        function plot_sd(obj,Vg)
            % draw the stability diagram
            
            N=zeros(1,length(Vg));           
            for j=1:length(Vg) 
                N(j)=obj.stable_conf(Vg(j));                    
            end;
            plot(Vg,N)
            
        end
        
        function Na=stable_conf(obj,V)
            % find stable charge configuration
            
            n1=0:obj.Nmax;
            E=zeros(1,obj.Nmax);
            
            for j1=1:length(n1)                   
                    E(j1)=obj.SD_SET([V;  0;   0], n1(j1));                   
            end;
            
            Na=n1(E(:)==min(E(:)));
            
        end
        
        function M=rate_matrix(obj,V)
            % computes a matrix of the transition rates for all possible
            % transitions which can occur in the system
            
            M=zeros(obj.Nmax, obj.Nmax);
            n=0:obj.Nmax-1;
            
            for j1=2:obj.Nmax-1
                for j2=1:obj.Nmax
                    
                    M(j1,j2)=obj.rates(n(j1),n(j2),'S', V)+obj.rates(n(j1),n(j2),'D', V);
                end;
            end;
                        
            M(1,2)=obj.rates(n(1),n(2),'S', V)+obj.rates(n(1),n(2),'D', V);
            M(obj.Nmax,obj.Nmax-1)=obj.rates(n(obj.Nmax),n(obj.Nmax-1),'S', V)+obj.rates(n(obj.Nmax),n(obj.Nmax-1),'D', V);
            
            for j=1:obj.Nmax
                M(j,j)=-sum(M(:,j));
            end;
            
        end
        
        function M=rate_matrix_s(obj,V)
            
            M=zeros(obj.Nmax, obj.Nmax);
            n=0:obj.Nmax-1;
            
            for j1=2:obj.Nmax-1
                for j2=1:obj.Nmax
                    
                    M(j1,j2)=obj.rates(n(j1),n(j2),'S', V);
                end;
            end;
                        
            M(1,2)=obj.rates(n(1),n(2),'S', V);
            M(obj.Nmax,obj.Nmax-1)=obj.rates(n(obj.Nmax),n(obj.Nmax-1),'S', V);
            
            for j=1:obj.Nmax
                M(j,j)=-sum(M(:,j));
            end;
            
        end        
        
    end
end
