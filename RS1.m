function dy = RS1(t,y,qd, vs)


%Vs=0.0215*(0.5*(1+tanh(4e10*(t-2e-9)))-0.5*(1+tanh(4e10*(t-4.0e-9))));        


dy=qd.rate_matrix([0.00 vs])*y;