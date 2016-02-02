function [Vg,Vs,Vd]=gate_volt(t)

omega=0.15e14;
A=0.003;
% Vg=0.0033+0.5*A*(1.0-cos(omega*t));
% %V=0.0033;




Vg=0.00831;


per1=18e-10;
start=5e-10;
per=[0 per1 2*per1 3*per1]+start;
%per1=9e-10;
%per=[0 per1 2*per1 3*per1 4*per1 5*per1 6*per1]+start;

for j=1:length(per)
    
    %Vg=Vg-A*(0.5*(1+tanh(4e10*(t-4e-10-per(j))))-0.5*(1+tanh(4e10*(t-6.5e-10-per(j)))))+0.00125*(0.5*(1+tanh(4e10*(t-4e-10-per(j))))-0.5*(1+tanh(4e10*(t-5.5e-10-per(j)))));
    Vg=Vg-A*(0.5*(1+tanh(4e10*(t-per(j))))-0.5*(1+tanh(4e10*(t-7e-10-per(j)))));
%     Vd=Vd-0.002*(0.5*(1+tanh(4e10*(t-4e-10-per(j))))-0.5*(1+tanh(4e10*(t-5.5e-10-per(j)))));        
%     Vd=Vd-0.002*(0.5*(1+tanh(4e10*(t-7.5e-10-per(j))))-0.5*(1+tanh(4e10*(t-9.0e-10-per(j)))));        
    %Vg=Vg+0.00325*0.7*(0.5*(1+tanh(4e10*(t-4e-10-per(j))))-0.5*(1+tanh(4e10*(t-6.5e-10-per(j)))));
    
end;

% per1=7e-10;
% per1=3.5e-10;
% per=[0 per1 2*per1 3*per1];
% per=[0 per1 2*per1 3*per1 4*per1 5*per1 6*per1];
% 
% Vg=0.0085;
% Vs=0.0;
% 
% for j=1:length(per)
%     
%     Vs=Vs+0.00025*(0.5*(1+tanh(4e10*(t-4e-10-per(j))))-0.5*(1+tanh(4e10*(t-6.0e-10-per(j)))));        
%     
% end;

% %Vg=0.0085-A*(0.5*(1+tanh(4e10*(t-4e-10)))-0.5*(1+tanh(4e10*(t-6.5e-10))))-A*(0.5*(1+tanh(4e10*(t-7.5e-10)))-0.5*(1+tanh(4e10*(t-10e-10))));
% %Vg=0.0085-A*(0.5*(1+tanh(4e14*(t-4e-10)))-0.5*(1+tanh(4e14*(t-8e-10))));
% %Vg=0.01-A*(0.5*(1+tanh(4e14*(t-4e-10)))-0.5*(1+tanh(4e14*(t-8e-10))));
% Vg=0;
% Vs=0;

%Vg=(0.0085+0*0.5*A)-A*(0.5*(1+tanh(5e10*(t-10e-10)))-0.5*(1+tanh(5e10*(t-30e-10))));
% Vd=0.0005*(0.5*(1+tanh(4e10*(t-4e-10)))-0.5*(1+tanh(4e10*(t-6e-10))));


Vs=0.0+t-t;
Vd=0.0+t-t;
%Vs=0.0+t-t;

%Vg=0.0049+0.0000001e14*(t-50e-14);
%Vs=0.0049+0.0000001e14*t;
