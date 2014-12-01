function g=ff(E,T)
    kB=1.38e-23;
    g=-E./(1-exp(E./kB./T));
end