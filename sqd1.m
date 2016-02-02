classdef sqd1
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
        
        function obj = sqd1(C,Cdv,Nmax,T)
            obj.RS=1000000;
            obj.RD=1000000;
                  obj.RS=1000000;
            obj.RD=1000000;      
%             obj.RS=4000000;
%             obj.RD=4000000;
            obj.num_states=2;
            % class constructor
            if (nargin == 0)
                q=1.6e-19;
                obj.T=0.03;
                obj.num_vn=3;
                obj.num_dots=1;
                obj.Nmax=15;
                obj.Cv=-[300  200  200].*q;                %cooupling between dots and voltage nodes
                obj.C=100*q-sum(obj.Cv);                              %cooupling between dots                
            elseif (nargin==2)
                obj.T=0.03;
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
            Vd=zeros(obj.num_dots,1);
            
            E=0.5*(obj.C\(-q.*N-obj.Cv*V))'*(-q.*N-obj.Cv*V); %free energy
            
            for j=1:obj.num_dots
                Naux=N(j)-1;
                mu_minus(j)=E-0.5*(obj.C\(-q.*Naux-obj.Cv*V))'*(-q.*Naux-obj.Cv*V);
                Naux=N(j)+1;
                mu_plus(j)=   0.5*(obj.C\(-q.*Naux-obj.Cv*V))'*(-q.*Naux-obj.Cv*V)-E;            
                
            end;
            
            E=E./q;
            Vd=obj.C\(-q.*N-obj.Cv*V);

            varargout{1} = mu_minus;
            varargout{2} = mu_plus;
            varargout{3} = Vd;
            
        end
        
%         function G=rates1(obj,i,j,V,i1,j1)
%             % the function computes the transition rate from charge state i
%             % to the charge state j at given voltages V
%             % flag specifies what lead is involved in the transition
%             
%             q=1.6e-19;            
%                                     
%             [~,mu_minus,mu_plus]=obj.SD_SET(V, i);
%             
%             if ((i-j)==1) %release an electron
%                 E=q*Vsd+mu_minus;
%             elseif ((i-j)==-1) %take an electron
%                 E=-q*Vsd-mu_plus;
%             else
%                 E=q;
%             end;
%             
%             if nargin>5 % there are several well defined energy states
%                 G(j1,j2)=1/q/R*ff_q(E,obj.T);
%             else        % there is several well defined energy states
%                 G=1/q/R/q*ff(E,obj.T);
%             end;
%             
%         end        
        
        
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
                E=-q*Vsd-mu_plus;
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
            switch obj.num_dots
                case 1
                    n1=n1;
                case 2
                    n1=combvec(n1,n1);
                case 3
                    n1=combvec(n1,n1,n1);
                case 4
                    n1=combvec(n1,n1,n1,n1);
                case 5
                    n1=combvec(n1,n1,n1,n1,n1);
                case 6
                    n1=combvec(n1,n1,n1,n1,n1,n1);
                case 7
                    n1=combvec(n1,n1,n1,n1,n1,n1,n1);   
                otherwise
                Disp('too many dots')        
            end;
            
            s=size(n1);
            E=zeros(1,s(2));
            
            for j1=1:s(2)                   
                    E(j1)=obj.SD_SET(V, n1(:,j1));                   
            end;
            
            Na=n1(:,E(:)==min(E(:)));
            Na=Na(:,1);
            
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
        
        function [I,Gamma1]=current(obj,vs,p)
            %             a=obj.states(1:end-1)-obj.states(2:end);
            %             b=find(a);
            %             num_ch_states=length(b)+1;
            q=1.6e-19;
            n=0:obj.Nmax-1;            
            I=zeros(1,length(squeeze(vs(1,:))));
            Gamma1=zeros(length(squeeze(vs(1,:))),obj.Nmax);
            
            for i=1:length(squeeze(vs(1,:)))
                
                Gamma=zeros(1,obj.Nmax);
                
                for j=1:obj.Nmax
                    for j1=1:obj.Nmax
                        
                        if j1-j==1
                            Gamma(j)=Gamma(j)+obj.rates(n(j1),n(j),'D', vs(:,i));
                        elseif j1-j==-1
                            Gamma(j)=Gamma(j)-obj.rates(n(j1),n(j),'D', vs(:,i));
                        else
                            
                        end;
                    end;
                end;          
                I(i)=q*sum(p(i,:).*Gamma);
                Gamma1(i,:)=Gamma;
            end;
        end
        
        function [I,p,Gamma] = current_st(obj,vs)
            %             a=obj.states(1:end-1)-obj.states(2:end);
            %             b=find(a);
            %             num_ch_states=length(b)+1;
            q=1.6e-19;
            n=0:obj.Nmax-1;
            I=zeros(1,length(squeeze(vs(1,:))));
            Gamma1=zeros(length(squeeze(vs(1,:))),obj.Nmax);
            
            for i=1:length(squeeze(vs(1,:)))
                Gamma=zeros(1,obj.Nmax);
                for j=1:obj.Nmax
                    for j1=1:obj.Nmax
                        if j1-j==1
                            Gamma(j)=Gamma(j)-obj.rates(n(j1),n(j),'D', vs(:,i));
                        elseif j1-j==-1
                            Gamma(j)=Gamma(j)+obj.rates(n(j1),n(j),'D', vs(:,i));
                        else
                            
                        end;
                    end;
                end;
                
                M=obj.rate_matrix(vs(:,i));
                M(1,:)=ones(1,obj.Nmax);
                B=zeros(obj.Nmax,1);
                B(1,1)=1;
                p=M\B;
                p=p';
                
                I(i)=q*sum(p.*Gamma);
                Gamma1(i,:)=Gamma;
                
            end;
            
        end           
        
    end
end
