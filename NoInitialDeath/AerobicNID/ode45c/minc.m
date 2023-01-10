function y=minc(z1,z2)
 if nargin <2
     ind=find(min(real(z1)));
     y=z1(ind);
     %y=min(z1,z2);
 else
        if (real(z1)<=real(z2));
            y=z1;
        else
            y=z2;
        end
        %y=min(z1,z2);
 end
end