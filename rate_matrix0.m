function G=rate_matrix0(VV,qd)
%  clear all
%  q=1.6e-19;
%  
%  
%  c1=300.*q;
%  c1s=500.*q;
%  
%  Cdv=-[c1s c1];
%  
%  C=-sum(Cdv);
%  qd=sqd1(C,Cdv);
% % 
% % % vector of states
% % 
% % V=-10:0.05:40;
% % Na=zeros(1,length(V));
% % 
% % for j=1:length(V)
% %     VV=[0*0.007e-3;  V(j)*1e-3];
% %     Na(j)=qd.stable_conf(VV);
% % end;
% 
% VV=[0; 0.021];
%--------------------------------------------------------

q=1.6e-19;
vs=0:20;
s=size(vs);
E=zeros(s(2),s(2));
G=zeros(s(2),s(2));
E1=zeros(s(2));

% for j=1:length(V)    
%              VV=[0;  V(j)*1e-3];
    for j1=1:s(2)

        E1(j1)=qd.SD_SET(VV, vs(j1));
        for j2=1:s(2)
%             E(j1,j2)=E1(j1)-qd.SD_SET(VV, vs(j2));
            if ((vs(j1)-vs(j2))==1) %release an electron
                E(j1,j2)=E1(j1)-qd.SD_SET(VV, vs(j2));
                G(j1,j2)=1/q/qd.RS/q*ff((q*(-0*VV(2)+E(j1,j2))),qd.T);
            elseif ((vs(j1)-vs(j2))==-1) %take an electron
                E(j1,j2)=-E1(j1)+qd.SD_SET(VV, vs(j2));
                G(j1,j2)=1/q/qd.RS/q*ff((q*(-0*VV(2)-E(j1,j2))),qd.T);
            else
                G(j1,j2)=0;
            end;
                

        end;
    end;
% end;

for j=1:s(2)    
    G(j,j)=0;
    G(j,j)=-sum(G(:,j));
end;

% plot(E1)