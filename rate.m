qd=sqd()

Vg=0.00:0.000025:0.02;

for j=1:length(Vg)
    G1(j)=qd.rates(2,3,'S', [0.007;  0.0;   0])+qd.rates(2,3,'D', [Vg(j);  0.0001;   0]);
    G2(j)=qd.rates(3,2,'S', [0.007;  0.0;   0])+qd.rates(3,2,'D', [Vg(j);  0.0001;   0]);
    G3(j)=qd.rates(2,1,'S', [0.007;  0.0;   0])+qd.rates(2,1,'D', [Vg(j);  0.0001;   0]);
    G4(j)=qd.rates(1,2,'S', [0.007;  0.0;   0])+qd.rates(1,2,'D', [Vg(j);  0.0001;   0]);
end;

figure(2)
plot(G1+G2+G3+G4)