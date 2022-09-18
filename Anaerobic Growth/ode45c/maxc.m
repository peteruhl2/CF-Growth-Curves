function y=maxc(z1,z2)
 if nargin <2
     ind=find(max(real(z1)));
     y=z1(ind);
 else
        if (real(z1)>=real(z2));
            y=z1;
        else
            y=z2;
        end
 end
end