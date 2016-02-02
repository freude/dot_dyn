function G=rate_matrix(VV,qd)

vs=combvec(1:7,1:7);
s=size(vs);
E=zeros(s(2),s(2));
G=zeros(s(2),s(2));

for j1=1:s(2)
    E1(j1)=qd.SD_SET(VV, vs(:,j1));
    for j2=1:s(2)
        E(j1,j2)=E1(j1)-qd.SD_SET(VV, vs(:,j2));
        if ((vs(2,j1)-vs(2,j2))~=0)&&((vs(1,j1)-vs(1,j2))==0)
            if ((vs(2,j1)-vs(2,j2))==1) %release an electron
                G(j1,j2)=1/q/qd.RS*ff_q(-(q*(VV(3)+E(j1,j2))),qd.T);
            elseif ((vs(2,j1)-vs(2,j2))==-1) %take an electron
                G(j1,j2)=1/q/qd.RS*ff_q((q*(VV(3)+E(j1,j2))),qd.T);
            else
                G(j1,j2)=0;
            end;
            G(j1,j2)=1/q/qd.RS*ff_q((q*(-VV(1)+E(j1,j2))),qd.T);
        elseif ((vs(2,j1)-vs(2,j2))==0)&&((vs(1,j1)-vs(1,j2))~=0)
            %             if ((vs(1,j1)-vs(1,j2))==1) %release an electron
            %                 G(j1,j2)=1/q/qd.RS*ff_q(-(q*(E(j1,j2))),qd.T);
            %             elseif ((vs(1,j1)-vs(1,j2))==-1) %take an electron
            %                 G(j1,j2)=1/q/qd.RS*ff_q((q*(E(j1,j2))),qd.T);
            %             else
            %                 G(j1,j2)=0;
            %             end;
            
            G(j1,j2)=1/q/qd.RS*ff_q((q*((E(j1,j2)))),qd.T);
        else
            G(j1,j2)=0;
        end;
    end;
end;

for j=1:s(2)    
    G(j,j)=0;
    G(j,j)=-sum(G(:,j));
end;