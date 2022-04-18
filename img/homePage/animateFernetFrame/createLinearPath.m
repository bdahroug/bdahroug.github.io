function L = createLinearPath(Pfirst, Psecond, iter)
    
    if (nargin<3)
      iter = 10;
    end
    
    deltaS = 0:1/iter:1;
    for i=1:iter+1
%         L(:,i) = (Psecond + l*deltaS*(Pfirst - Psecond))';
        L(:,i) = Pfirst + deltaS(i)*(Psecond - Pfirst);
    end    
    
end