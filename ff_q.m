function g=ff_q(E,T)
    kB=1.38e-23;
    g=1/(1+exp(E./kB./T));
end