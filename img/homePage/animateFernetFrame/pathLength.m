%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute curvilinear length (path length) _ version 1
% created by Jean-Antoine SEON 
% edited by Bassem DAHROUG _ 29-01-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I/P: 
%     - curve(:,3) = [xp; yp; zp] 
%     
%  O/P:
%     - curve(:,4) = [xp; yp; zp; s] , s: path length
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = pathLength(M)
	
	l = length(M) ;
% 	[m n] = size(M);
	M(4,1) = 0 ;     
	
    for i=2:l
%         deltaS = norm( M(1:3,i) - M(1:3,i-1) )
        M(4,i) = M(4,i-1) + norm( M(1:3,i) - M(1:3,i-1) );
    end
    
    
	