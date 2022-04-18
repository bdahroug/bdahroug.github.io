%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute arc curvature _ version 1
% edited by Bassem DAHROUG _ 29-01-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I/P: 
%     - curve(:,4) = [xp; yp; zp; s] 
%     
%  O/P:
%     - curve(:,5) = [xp; yp; zp; s; c] , c:path curvature
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = pathCurvature(M)

	l = length(M) ;
% 	[m n] = size(M);
%     
%     %first intersection point 
%     normalVector = cross( (M(1:3,1) - M(1:3,l)), (M(1:3,2) - M(1:3,1)) ); %normal vector on the plan 
%     n = normalVector/ norm ( normalVector );  
%     Pm = ( M(1:3,1) + M(1:3,l) ) / 2; % middle point    
%     a = [n'; n'; (M(1:3,2) - M(1:3,1))'] ;
%     b = [n'*M(1:3,1); n'*M(1:3,2); (M(1:3,2) - M(1:3,1))'*Pm];
%     P_inter = a \ b;        % intersection point 
%     %curvarure radius
%     checkA=
%     if (norm(b) == 0) || (~isnan(a))
%         M(5,1) = 0;
%     else
%         M(5,1) = 1/norm( Pm - P_inter(1:3) );
%     end
%     
    M(5,1) = 0;
    for i=2:l
        if i==l-1
            v1 = M(1:3,i+1) - M(1:3,i);
            v2 = M(1:3,1) - M(1:3,i);
        elseif i==l
            v1 = M(1:3,1) - M(1:3,i);
            v2 = M(1:3,2) - M(1:3,i);
        else
            v1 = M(1:3,i+1) - M(1:3,i);
            v2 = M(1:3,i+2) - M(1:3,i);
        end
        n = cross(v1, v2)/ norm ( cross(v1, v2) )  %normal vector on the plan 
        
        if isnan(norm(n))%straight line case
            M(5,i) = 0;
        else
            Pm = ( M(1:3,i) + M(1:3,i-1) ) / 2; % middle point    
            a = [n'; n'; v1'] ;
            if i==l
                b = [n'*M(1:3,i); n'*M(1:3,1); v1'*Pm];
            else
                b = [n'*M(1:3,i); n'*M(1:3,i+1); v1'*Pm];
            end
            P_inter = a \ b;        % intersection point 
            %curvature radius
            if isnan(norm(P_inter)) %isnan(norm(b))
                M(5,i) = 0;
            else
                M(5,i) = 1/norm( Pm - P_inter(1:3) );
            end
        end
    end
    
%
%     %first intersection point (u^v).w=0, u=P1P2, v=P1Pi, w=PmPi
%     L12 = M(1:3,1) - M(1:3,l);   %vector line between the first & second point
%     skew_L12 = skew(L12);
%     Pm = ( M(1:3,1) + M(1:3,l) ) / 2; % middle point
%     b= cross(L12,M(1:3,l)) + Pm;
%     P_inter = skew_L12\b(1:3);        % intersection point 
%     %curvarure radius
%     if det(skew_L12) == 0
%         M(5,1) = 0;
%     else
%         M(5,1) = 1/norm( Pm - P_inter(1:3) );
%     end
%     
%     for i=2:l
%         L = M(1:3,i) - M(1:3,i-1);   %vector line between the first & second point
%         skew_L = skew(L);
%         Pm = (M(1:3,i) + M(1:3,i-1)) / 2; % middle point
%         b= cross(L,M(1:3,i-1)) + Pm;
%         P_inter = skew_L\b;        % intersection point 
%         %curvarure radius
%         if det(skew_L) == 0
%             M(5,i) = 0;
%         else
%         	M(5,i) = 1/norm( M(1:3,i) - P_inter(1:3) );
%         end
%     end
    

end 