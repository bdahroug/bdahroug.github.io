%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute a geometric curve in 3D 
% created by Bassem DAHROUG _ 06-01-17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I/P: 
%     - starting point of curve Ps=[xs; ys; zs]
%     - penetration point Pp=[xp; yp; zp]
%     - unitType = 0(mm) & 1(m)
%     - curveType = 0(helical) & 1(zigzag) & 2(uros)
%     
%  O/P:
%     - curve(5,:) = [xp; yp; zp; s] , s:path length 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function curve = geometricCurve(pathStartPoint, penetrationPoint, unitType, curveType)
    
if unitType==1
    firstLineLength=0.005;
else
    firstLineLength=10;
end

iter1=10;
L1 = createLinearPath(penetrationPoint-[0;firstLineLength;0], penetrationPoint, iter1);
L2 = createLinearPath(penetrationPoint, pathStartPoint, 2*iter1);
    

init__w_curve= [expm(skew([0;0;0]*pi/180)), pathStartPoint; 0 0 0 1];
pathStartPoint=zeros(3,1);
% figure, hold on;
% trace_repere(init__w_curve);

w_R_curveShifted = expm(skew([0;0;0]*pi/180));
% w_R_curveShifted = expm(skew([0;0;45]*pi/180));
w_P_curveShifted = [0;0;0];
w_T_curveShifted = init__w_curve*[w_R_curveShifted, w_P_curveShifted; 0 0 0 1];
% trace_repere(w_T_curveShifted);


%generated curve
if curveType==0     %spiral curve
    iter2=150;
    if unitType==1  %m
        r_max=0.005; helixHeightGain=0.001;
    else    %mm
        r_max=20; helixHeightGain=1;
    end
    r=0:r_max/iter2:r_max;
    t=0:6*pi/iter2:6*pi;
    x = r.*cos(t) + pathStartPoint(1);
    y = helixHeightGain*t + pathStartPoint(2);
    z = r.*sin(t) + pathStartPoint(3);
    curveSpiral = [x; y; z];
%     curve = [L1 L2(:,2:end) curveSpiral(:,2:end)];
    w_curvePointsShifted = w_T_curveShifted*[curveSpiral; ones(1,length(curveSpiral))];  
    curve = [L1 L2(:,2:end) w_curvePointsShifted(1:3,2:end)];
    
    
elseif curveType==1        %zigzag curve
    iter2=30; iter3=15;
    point1 = pathStartPoint+[-0.006;0.008;0.003];
    point2 = pathStartPoint+[0;0.01;0.005];
    if unitType==0 %mm
        point1 = point1*100;
        point2 = point2*100;
    end
    L3= createLinearPath(pathStartPoint, point1, iter2);
    L4= createLinearPath(point1, point2, iter3);
    curveZigzag = [L3(:,2:end) L4(:,2:end)];
%     curve = [L1 L2(:,2:end) curveZigzag(:,2:end)];
    w_curvePointsShifted = w_T_curveShifted*[curveZigzag; ones(1,length(curveZigzag))];  
    curve = [L1 L2(:,2:end) w_curvePointsShifted(1:3,2:end)];
    
    
elseif curveType==2     %uRoCs curve
    iter2=20;
    u=[0,  60, 60, 120, 120, 200, 200, 160, 200, 240, 320, 280, 240, 320, 400, 360, 360, 400, 540, 540, 460, 460, 540];
    v=[0, 160, 40,  40, 160, 160, 120, 120,  40,  40, 100, 160, 100,  40,  40,  80, 120, 160, 160, 100, 100,  40,  40];
    w=zeros(1,length(u));
    curvePoints= [u;v;w];
    if unitType==1
        curvePoints= curvePoints./1000/10;
    else
        curvePoints= curvePoints./10;
    end

    angleDeg=[0;-45;0];
    w_R_curve= expm(skew(angleDeg*pi/180));
    w_P_curve= [0;0;0];
    w_T_curve = [w_R_curve,w_P_curve; 0 0 0 1];
    w_curvePoints = w_T_curve*[curvePoints; ones(1,length(u))];

    curveUros=[];
    for i=2:length(u)
        L_uros= createLinearPath(pathStartPoint+w_curvePoints(1:3,i-1), pathStartPoint+w_curvePoints(1:3,i), iter2);
        curveUros= [curveUros L_uros(:,2:end)];
    end
%     curve = [L1 L2(:,2:end) curveUros(:,2:end)];
    w_curvePointsShifted = w_T_curveShifted*[curveUros; ones(1,length(curveUros))];  
    curve = [L1 L2(:,2:end) w_curvePointsShifted(1:3,2:end)];
    
elseif curveType==3     %curved path

    circleIter=10; 
    if unitType==1  %m
        circleRadiusX=0.0025; circleRadiusY=0.01;
    else    %mm
        circleRadiusX=25; circleRadiusY=10;
    end
    lineLength= norm(pathStartPoint-penetrationPoint);
    circleCentre(:,1)= penetrationPoint + [-circleRadiusX;lineLength+circleRadiusY;0];
    
    thetaMin=0; thetaMax=120*pi/180; deltaTheta= (thetaMax-thetaMin)/circleIter;
    theta= thetaMin:deltaTheta:thetaMax;
    x_t_circle= circleRadiusX*cos(theta) + circleCentre(1);
    y_t_circle= circleRadiusY*sin(theta) + circleCentre(2);
    z_t_circle= circleCentre(3)*ones(1, length(x_t_circle));
    
    circle1=[x_t_circle; y_t_circle; z_t_circle]; 
    curve = [L1 L2(:,2:end) circle1(:,2:end)];
    
end
    

curve = pathLength(curve);
%curve = pathCurvature(curve);

% % figure, hold on;
% % scatter(curve(:,1), curve(:,2), '+');
% plot3(curve(1,:), curve(2,:), curve(3,:),'.');
% xlabel('x'); ylabel('y'); zlabel('z');
% grid on, hold off;

end

