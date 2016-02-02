clear all

qd=sqd1()

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

[T,p] = ode45(@(t,y) RS(t,y,qd), [0 10.0e-9], [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0]);
% 
[V1,V2,V3]=gate_volt(T);

figure(1)
hold on;
%plot(T,p(:,:),T,V1./max(abs(V1)),T,V2./max(abs(V2)),T,V3./max(abs(V3)))
plot(T,p(:,:))
figure(3)
hold on
plot(T,V1,T,V2,T,V3)
% 
%[I,Gamma]=qd.current([V1';V2';V3'],p);
%[I,p1,Gamma]=qd.current_st([V1';V2';V3']);

Vg=0.005:0.0001:0.012;
Vsd=-0.0015:0.00005:0.0015;

for j1=1:length(Vg)
    j1
    for j2=1:length(Vsd)
        [I(j1,j2),p1(j1,j2,:),Gamma(j1,j2,:)]=qd.current_st([Vg(j1); Vsd(j2); 0]);
        r(j1,j2)=qd.rates(3,2,'D', [Vg(j1); Vsd(j2); 0])+qd.rates(3,2,'S', [Vg(j1); Vsd(j2); 0]);
        r1(j1,j2)=qd.rates(2,3,'D', [Vg(j1); Vsd(j2); 0])+qd.rates(2,3,'S', [Vg(j1); Vsd(j2); 0]);
    end;
end;

%     
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