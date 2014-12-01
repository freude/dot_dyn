clear all

qd=rate_model()

% Vg=-0.02:0.000025:0.03;
% 
% for j=1:length(Vg)
%     G1(j)=qd.rates(1,3,'S', [Vg(j);  0.0]);%+qd.rates(1,2,'D', [Vg(j);  0.0]);
%     G2(j)=qd.rates(3,1,'S', [Vg(j);  0.0]);%+qd.rates(2,1,'D', [Vg(j);  0.0]);
% %     G3(j)=qd.rates(4,6,'S', [Vg(j);  0.0]);%+qd.rates(1,2,'D', [Vg(j);  0.0]);
% %     G4(j)=qd.rates(6,4,'S', [Vg(j);  0.0]);%+qd.rates(2,1,'D', [Vg(j);  0.0]);
% end;
% 
% figure(2)
% plot(Vg,G1,Vg,G2)

%[T,p] = ode45(@(t,y) RS(t,y,qd), [0 1150e-12], [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0]);

% vs=-1.5;
% 
% [T,p] = ode45(@(t,y) RS1(t,y,qd,vs), [0 2e-9], [0 1 0 0]);
% 
% hold on
% plot(T,p(:,:))


vs=0:-0.01:-1.5;
I=qd.current_st(vs);
hold on
plot(vs(1:end-1),diff(I)./diff(vs));

% for j=1:length(vs)
%     rate1(j)=qd.rates(2,1,'S', [0.0 vs(j)]);
%     rate2(j)=qd.rates(3,1,'S', [0.0 vs(j)]);
%     rate3(j)=qd.rates(4,1,'S', [0.0 vs(j)]);
% end;
% 
% hold on
% plot(vs,rate1,vs,rate2,vs,rate3)


% I=zeros(1,length(T));
% 
% for j=1:length(T)
%     M=qd.rate_matrix_s([V1(j);V2(j);V3(j)]);
%     for j1=2:qd.Nmax-1
%         I(j)=I(j)+p(j,j1)*(M(j1+1,j1)-M(j1-1,j1));
%     end;
% end;
% 
% figure(2)
% plot(T,p(:,:),T,V1./max(V1),T,V2./max(V2))

