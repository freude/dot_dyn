clear all

qd=sqd() 

Vg=0:0.0003:0.07;
 
 
 for j=1:length(Vg)
     r1(j)=qd.rates(2,3,'S',[Vg(j);0.001;0])
     r2(j)=qd.rates(3,2,'S',[Vg(j);0.001;0])
 end;
 
 figure(2)
 plot(Vg,r1,Vg,r2)