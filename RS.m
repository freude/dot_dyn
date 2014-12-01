function dy = RS(t,y,qd)


[V1,V2,V3]=gate_volt(t);

dy=qd.rate_matrix([V1;V2;V3])*y;

