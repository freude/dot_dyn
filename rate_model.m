classdef rate_model
    % the class implements the quantum dot system
    
    properties
        % define the properties of the class here, (like fields of a struct)
        spectr;
        states;
        Nmax;
        T;
        J_tip;
        J_s;
        x;
        L
        alpha;
    end
    
    methods
        % methods, including the constructor are defined in this block
        
        function obj = rate_model()            
            q=1.6e-19;            
            obj.J_tip=semicond_metal();
            obj.J_s=semicond_semicond();
            obj.T=4;            
            obj.spectr=[0 5 9.2 11.2].*q.*0.001;   % energy spectrum relative to the Fermi energy
            obj.states=[1 2 2 2];   %charge state vs excitation spectrum    
            obj.x=5e-9;
            obj.L=10e-9;
            obj.alpha=0.1;            
        end
               
        function G=rates(obj,j1,j2,flag,V)
            % the function computes the transition rate from charge state i
            % to the charge state j at given voltages V
            % flag specifies what lead is involved in the transition
            
            q=1.6e-19;            
            
            vv=(V(2)-V(1))*obj.alpha*(1-obj.x/obj.L);

            
            if obj.states(j1)==obj.states(j2)
                dE=j2-j1;
                if (dE==1)&&(j1==2)
                    G=0*1e13;
                elseif (dE==1)&&(j1==3)
                    G=0*2e13;                    
                elseif (dE==2)&&(j1==2)
                    G=0*1e13;                    
                else
                    G=0;
                end;
                return
                
            end;            
            
            if flag=='S'
                
                Vsd=V(1);
                d=5e-9;
                
                if obj.states(j1)>obj.states(j2)                    
                    dE=obj.spectr(j1)-obj.spectr(j2)-q*Vsd+q*vv;
                    G=0.1e14*obj.J_s.transparency(d,obj.spectr(j1),0)*ff_q(dE,obj.T);                    
                else                    
                    dE=-obj.spectr(j2)+obj.spectr(j1)+q*Vsd-q*vv;
                    G=0.1e14*obj.J_s.transparency(d,obj.spectr(j2),0)*ff_q(dE,obj.T);
                end;               
                
            elseif flag=='D'
                
                Vsd=V(2);
                d=0.01e-9;
                
                if obj.states(j1)>obj.states(j2) % give an electron                    
                    dE=obj.spectr(j1)-obj.spectr(j2)-q*Vsd+q*vv;
                    G=1e12*obj.J_tip.transparency(d,obj.spectr(j1)-obj.spectr(j2)+q*vv,Vsd)*ff_q(dE,obj.T);                    
                else                           % take an electron                    
                    dE=-obj.spectr(j2)+obj.spectr(j1)+q*Vsd-q*vv;
                    G=1e12*obj.J_tip.transparency(d,-obj.spectr(j2)+obj.spectr(j1)-q*vv,Vsd)*ff_q(dE,obj.T);                    
                end;                                
            end;            
            
        end
        
        function M=rate_matrix(obj,V)
            % computes a matrix of the transition rates for all possible
            % transitions which can occur in the system
            
            M=zeros(length(obj.spectr), length(obj.spectr));
            
            for j1=1:length(obj.spectr)
                for j2=1:length(obj.spectr)
                    if obj.states(j1)==obj.states(j2)
                        M(j1,j2)=obj.rates(j1,j2,'R', V);
                    else
                        M(j1,j2)=obj.rates(j1,j2,'S', V)+obj.rates(j1,j2,'D', V);
                    end;
                end;
            end;
            
            for j=1:length(obj.spectr)
                M(j, j)=-sum(M(:,j));
            end;
                        
            
        end
        
        function I=current(obj,vs)
            %             a=obj.states(1:end-1)-obj.states(2:end);
            %             b=find(a);
            %             num_ch_states=length(b)+1;
            q=1.6e-19;
            
            I=zeros(1,length(vs));
            
            for i=1:length(vs)
                
                Gamma=zeros(1,length(obj.spectr));
                
                for j=1:length(obj.spectr)
                    for j1=1:length(obj.spectr)
                        
                        if (obj.states(j1)-obj.states(j))==1
                            Gamma(j)=Gamma(j)+obj.rates(j1,j,'S', [0.0 vs(i)]);
                        elseif (obj.states(j1)-obj.states(j))==-1
                            Gamma(j)=Gamma(j)-obj.rates(j1,j,'S', [0.0 vs(i)]);
                        else
                            
                        end;
                    end;
                end;
                
                [~,p] = ode45(@(t,y) RS1(t,y,obj,vs(i)), [0 4.0e-9], [1 0 0 0]);
                I(i)=q*sum(p(end,:).*Gamma);
            end;            
        end
        
        function [I,p1,Gamma1] = current_st(obj,vs)
            %             a=obj.states(1:end-1)-obj.states(2:end);
            %             b=find(a);
            %             num_ch_states=length(b)+1;
            q=1.6e-19;
            
            I=zeros(1,length(vs));
            p1=zeros(4,length(vs));
            Gamma1=zeros(4,length(vs));
            
            for i=1:length(vs)
                
                Gamma=zeros(1,length(obj.spectr));
                
                for j=1:length(obj.spectr)
                    for j1=1:length(obj.spectr)                        
                        if (obj.states(j1)-obj.states(j))==1
                            Gamma(j)=Gamma(j)+obj.rates(j1,j,'D', [0.0 vs(i)]);
                        elseif (obj.states(j1)-obj.states(j))==-1
                            Gamma(j)=Gamma(j)+obj.rates(j1,j,'D', [0.0 vs(i)]);
                        else
                            
                        end;
                    end;
                end;
                
                M=obj.rate_matrix([0.00 vs(i)]);
                M(1,:)=ones(1,length(obj.spectr));
                B=zeros(length(obj.spectr),1);
                B(1,1)=1;
                p=M\B;
                p=p';
                
                I(i)=q*sum(p.*Gamma);
                p1(:,i)=p;
                Gamma1(:,i)=Gamma;
            end;            
        end        
        
        
    end
end
