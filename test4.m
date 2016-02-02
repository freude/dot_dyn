clear all
q=1.6e-19;

c10=300.*q;
c20=300.*q;
Cdv=[c10 0; 0 c20];

c12=50.*q;
c1s=50.*q;
c2s=50.*q;

C=[c1s+c12+Cdv(1,1) c12; c12 c2s+c12+Cdv(2,2)];


qd=sqd1(C,Cdv);



V=0:0.4:15;
m1=zeros(length(V),length(V),2);
m2=zeros(length(V),length(V),2);
Vd=zeros(length(V),length(V),2);
N=zeros(length(V),length(V),2);




for j1=1:length(V)
    for j2=1:length(V)
        VV=[-V(j1)*1e-3;  -0*V(j2)*1e-3];
    [E1 ,m11 , m21, Vd1] =qd.SD_SET(VV, [2; 2]);
    E(j1,j2)=E1;
    m1(j1,j2,:)=m11;
    m2(j1,j2,:)=m21;
    Vd(j1,j2,:)=Vd1;
    
    N1=qd.stable_conf(VV);
    N(j1,j2,1)=N1(1);
    N(j1,j2,2)=N1(2);
    
    end;
end;

surf(squeeze(N(:,:,1))+squeeze(N(:,:,2)))

