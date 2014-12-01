 function [Na,varargout]=SD_SET(V,C,Cv,Num,T)

kb=1.38e-23;
q=1.6e-19;
n1=0:Num;
dT=0.5*kb*T;
I=0;

E=zeros(length(n1),length(n1));
ch_pot1=zeros(length(n1),length(n1));
ch_pot2=zeros(length(n1),length(n1));
ch_pot3=zeros(length(n1),length(n1));
ch_pot4=zeros(length(n1),length(n1));

Na=[];

for j1=1:length(n1)
    for j2=1:length(n1)
                
                E(j1,j2)=0.5*(C\(-q.*[n1(j1); n1(j2)]-Cv*V))'*(-q.*[n1(j1); n1(j2)]-Cv*V); %free energy
                
                ch_pot1(j1,j2)=E(j1,j2)-0.5*(C\(-q.*[n1(j1)-1; n1(j2)]-Cv*V))'*(-q.*[n1(j1)-1; n1(j2)]-Cv*V);
                ch_pot2(j1,j2)=0.5*(C\(-q.*[n1(j1)+1; n1(j2)]-Cv*V))'*(-q.*[n1(j1)+1; n1(j2)]-Cv*V)-E(j1,j2);
                
                ch_pot3(j1,j2)=E(j1,j2)-0.5*(C\(-q.*[n1(j1); n1(j2)-1]-Cv*V))'*(-q.*[n1(j1); n1(j2)-1]-Cv*V);
                ch_pot4(j1,j2)=0.5*(C\(-q.*[n1(j1); n1(j2)+1]-Cv*V))'*(-q.*[n1(j1); n1(j2)+1]-Cv*V)-E(j1,j2);
                                        
    end;
end;

% minimization of free energy to find cahrge configuration

i=find(E(:)==min(E(:)));

% processing several minima

ij=sqrt(-1);

if (length(i)>1)
    [m1,m2] = ind2sub(size(E),i(1));
%     [m11,m21,m31] = ind2sub(size(E),i(2));
%     Na=[n1(m1)+ij*n1(m11); n1(m2)+ij*n1(m21); n1(m3)+ij*n1(m31)];
    Na=[n1(m1); n1(m2)];
else
    [m1,m2] = ind2sub(size(E),i);
    Na=[n1(m1); n1(m2)];
end;

%[Vd,i] = min(E(:));

% find stability regions:
% flag=0 the region is totaly stable
% flag=1 dot only is unstable
% flag=2 sensor is unstable only
% flag=3 both are unstable

flag=0;
if (V(3)<0)
    con=((ch_pot1(m1,m2)>(V(3)*q-dT))||(ch_pot2(m1,m2)<(0+dT)));
    con2=((ch_pot1(m1,m2)>(V(3)*q-dT))&&(ch_pot1(m1,m2)<(0+dT)));
    con3=((ch_pot2(m1,m2)>(V(3)*q-dT))&&(ch_pot2(m1,m2)<(0+dT)));
else
    con=((ch_pot1(m1,m2)>(0-dT))||(ch_pot2(m1,m2)<(V(3)*q+dT)));
    con2=((ch_pot1(m1,m2)>(0-dT))&&(ch_pot1(m1,m2)<(V(3)*q+dT)));
    con3=((ch_pot2(m1,m2)>(0-dT))&&(ch_pot2(m1,m2)<(V(3)*q+dT)));
end;
if (V(4)<0)
    con1=((ch_pot3(m1,m2)>(V(4)*q-dT))||(ch_pot4(m1,m2)<(0+dT)));
    con4=((ch_pot3(m1,m2)>(V(4)*q-dT))&&(ch_pot3(m1,m2)<(0+dT)));
    con5=((ch_pot4(m1,m2)>(V(4)*q-dT))&&(ch_pot4(m1,m2)<(0+dT)));
else
    con1=((ch_pot3(m1,m2)>(0-dT))||(ch_pot4(m1,m2)<(V(4)*q+dT)));    
    con4=((ch_pot3(m1,m2)>(0-dT))&&(ch_pot3(m1,m2)<(V(4)*q+dT)));
    con5=((ch_pot4(m1,m2)>(0-dT))&&(ch_pot4(m1,m2)<(V(4)*q+dT)));
    
end;



if con
    flag=1;
end;

if con1
    flag=flag+2;
end;


if (flag==1)||(flag==3)
     
    if con2
        Na(1)=Na(1)-round(rand());I=1;
    end;
    if con3
        Na(1)=Na(1)+round(rand());I=1;
    end;
end;

A=0;

if (flag==2)
    
    if con4        
        N2=Na(2)-1;
        A=backaction2(V,C,Cv,Num,Na,N2);
        if (A==-1)
            Na(1)=Na(1)-round(rand());I=1;
        end;
        if (A==1)
            Na(1)=Na(1)+round(rand());I=1;
        end;      
        %Na(2)=Na(2)-round(rand());
    end;
    
    if con5
        N2=Na(2)+1;
        A=backaction2(V,C,Cv,Num,Na,N2);
        if (A==-1)
            Na(1)=Na(1)-round(rand());I=1;
        end;
        if (A==1)
            Na(1)=Na(1)+round(rand());I=1;
        end;       
        %Na(2)=Na(2)+round(rand());        
    end;
end;

Na(Na<0)=0;

minE=0.5*(C\(-q.*Na-Cv*V))'*(-q.*Na-Cv*V);

Nam=[Na(1)-1; Na(2)];
ch_pot1=minE-0.5*(C\(-q.*Nam-Cv*V))'*(-q.*Nam-Cv*V);
Nam=[Na(1)+1; Na(2)];
ch_pot2=0.5*(C\(-q.*Nam-Cv*V))'*(-q.*Nam-Cv*V)-minE;
Nam=[Na(1); Na(2)-1];
ch_pot3=minE-0.5*(C\(-q.*Nam-Cv*V))'*(-q.*Nam-Cv*V);
Nam=[Na(1); Na(2)+1];
ch_pot4=0.5*(C\(-q.*Nam-Cv*V))'*(-q.*Nam-Cv*V)-minE;

Vd=C\(-q.*Na-Cv*V);

varargout{1} = Vd;
varargout{2} = ch_pot1;
varargout{3} = ch_pot2;
varargout{4} = ch_pot3;
varargout{5} = ch_pot4;
varargout{6} = I;
varargout{7} = flag;
varargout{8} = minE;
varargout{9} = A;

