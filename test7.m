clear all
q=1.6e-19;

c1=300.*q;
c1s=500.*q;


Cdv=-[c1s c1];

C=-sum(Cdv);
qd=sqd1(C,Cdv);

% vector of states

% V=0:0.4:30;
% V1=0.0e-3;
% VV=[0;  V1];


% s=size(vs);
% E=zeros(s(2),s(2));
% G=zeros(s(2),s(2));
% G1=zeros(s(2),s(2));
% G2=zeros(s(2),s(2));
% vs=combvec(1:7,1:7);
%
% for j1=1:s(2)
%     E1(j1)=qd.SD_SET(VV, vs(:,j1));
%     for j2=1:s(2)
%         E(j1,j2)=E1(j1)-qd.SD_SET(VV, vs(:,j2));
%         if ((vs(2,j1)-vs(2,j2))~=0)&&((vs(1,j1)-vs(1,j2))==0)
%             if ((vs(2,j1)-vs(2,j2))==1) %release an electron
%                 G(j1,j2)=1/q/qd.RS*ff_q(-(q*(VV(3)+E(j1,j2))),qd.T);
%             elseif ((vs(2,j1)-vs(2,j2))==-1) %take an electron
%                 G(j1,j2)=1/q/qd.RS*ff_q((q*(VV(3)+E(j1,j2))),qd.T);
%             else
%                 G(j1,j2)=0;
%             end;
%             G(j1,j2)=1/q/qd.RS*ff_q((q*(-VV(1)+E(j1,j2))),qd.T);
%         elseif ((vs(2,j1)-vs(2,j2))==0)&&((vs(1,j1)-vs(1,j2))~=0)
%             %             if ((vs(1,j1)-vs(1,j2))==1) %release an electron
%             %                 G(j1,j2)=1/q/qd.RS*ff_q(-(q*(E(j1,j2))),qd.T);
%             %             elseif ((vs(1,j1)-vs(1,j2))==-1) %take an electron
%             %                 G(j1,j2)=1/q/qd.RS*ff_q((q*(E(j1,j2))),qd.T);
%             %             else
%             %                 G(j1,j2)=0;
%             %             end;
%
%             G(j1,j2)=1/q/qd.RS*ff_q((q*((E(j1,j2)))),qd.T);
%         else
%             G(j1,j2)=0;
%         end;
%     end;
% end;
%
% for j=1:s(2)
%     G(j,j)=0;
%     G(j,j)=-sum(G(:,j));
% end;

%figure(2);surf(G1)
% Vg=0:0.000025:0.02;

% for j=1:length(Vg)
%     G7(j)=qd.rates(4,3,'S', [Vg(j);  0.0;   0])+qd.rates(4,3,'D', [Vg(j);  0.01;   0]);
%     G8(j)=qd.rates(3,4,'S', [Vg(j);  0.0;   0])+qd.rates(3,4,'D', [Vg(j);  0.01;   0]);
%     G1(j)=qd.rates(2,3,'S', [Vg(j);  0.0;   0])+qd.rates(2,3,'D', [Vg(j);  0.01;   0]);
%     G2(j)=qd.rates(3,2,'S', [Vg(j);  0.0;   0])+qd.rates(3,2,'D', [Vg(j);  0.01;   0]);
%     G3(j)=qd.rates(2,1,'S', [Vg(j);  0.0;   0])+qd.rates(2,1,'D', [Vg(j);  0.01;   0]);
%     G4(j)=qd.rates(1,2,'S', [Vg(j);  0.0;   0])+qd.rates(1,2,'D', [Vg(j);  0.01;   0]);
%     G5(j)=qd.rates(0,1,'S', [Vg(j);  0.0;   0])+qd.rates(0,1,'D', [Vg(j);  0.01;   0]);
%     G6(j)=qd.rates(1,0,'S', [Vg(j);  0.0;   0])+qd.rates(1,0,'D', [Vg(j);  0.01;   0]);
% end;
%
% figure(2)
% plot(Vg,G2)
%
%[T,p] = ode45(@(t,y) RS(t,y,qd), [0 1150e-12], [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);

ic=zeros(1,21);
ic(1,7)=1;

%V1=-5e-3:1e-4:5e-3;
[T,p] = ode45(@(t,y) RS(t,y,qd,[0; 0]), [0 5.0e-9], ic);

hold on
plot(T,p)