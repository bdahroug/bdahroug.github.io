
clc, clear all, close all;

%% I/P
metricUnitType = 0;    %0(mm) & 1(m) 
curveType = 0;         %0(spiral) & 1(zigzag) & 2(uros) & 3(curved) & 4(lissajous)

% O/P
op_animatedGIF = 0;
op_saveImagesPNG = 0;
op_animatedAVI = 0;

%% generate the curve
init_w_P_r = [0; 0; 0; 1];    

if curveType<4
    if metricUnitType==0
        pathStartpoint = init_w_P_r(1:3) + [0; 10; 0];
    end
    w_refPathCurve = geometricCurve(pathStartpoint, init_w_P_r(1:3), metricUnitType, curveType);   

elseif curveType==4    % lissajous curve
    maxCurveSamplePoints = 100;
    curveWidth = 1; 
    curveHeight = 2;
    horizonLobe = 1; 
    verticalLobe = 2;
    deltaRotationAngle = ( (curveHeight-1)/curveHeight ) * (pi/2);
    stepsIter = linspace(0,10,maxCurveSamplePoints);
    scaleFactor = 4;

    x_lissajous = scaleFactor*curveWidth*sin((horizonLobe*stepsIter)+deltaRotationAngle);
    y_lissajous = scaleFactor*curveHeight*cos(verticalLobe*stepsIter);
    z_lissajous = zeros(maxCurveSamplePoints);

    %curveLissajous= [x_lissajous(1:62)', y_lissajous(1:62)', z_lissajous(1:62)', ones(1,62)']
    curveLissajous = [x_lissajous(24:62)', y_lissajous(24:62)', z_lissajous(24:62)', ones(1,39)';...
                      x_lissajous(1:24)', y_lissajous(1:24)', z_lissajous(1:24)', ones(1,24)'];

    w_refPathCurve = curveLissajous';
end

%% compute Fernet frame
epsilon = 1e-7;
lengthPath = length(w_refPathCurve);
T = zeros(3, lengthPath);       % unit tangent vector
N = zeros(3, lengthPath);       % principal normal vector
B = zeros(3, lengthPath);       % binormal vector
k = zeros(3, lengthPath);       % curvature
t = zeros(1, lengthPath);       % torsion of curve

for i=1:lengthPath
    if i==1
        %T(:,i) = (w_refPathCurve(1:3,i) - zeros(3,1)) / (w_refPathCurve(4,i) - 0);
    else
        delta__T(:,i) = w_refPathCurve(1:3,i) - w_refPathCurve(1:3,i-1);
        %T(:,i) = ( delta__T(:,i) ) / ( w_refPathCurve(4,i) - w_refPathCurve(4,i-1) );
        delta__norm_T(1,i) = norm( delta__T(:,i) );
        T(:,i) = delta__T(:,i) / delta__norm_T(1,i);
        
        %delta__norm_T(1,i) = norm(T(1:3,i)) - norm(T(1:3,i-1));
        delta__T_minus(:,i) = T(1:3,i) - T(1:3,i-1);
        delta__norm_T_minus(1,i) = norm( delta__T_minus(:,i) );
        if delta__norm_T_minus(1,i) < epsilon
            %N(:,i) = zeros(3, 1);
            N(:,i) = [-1;0;0];
        else
            N(:,i) = delta__T_minus(:,i) / ( delta__norm_T_minus(1,i) );
        end
        
        TxN(:,i) = cross( T(:,i), N(:,i) );
        norm_TxN(1,i) = norm(TxN(:,i));
        if norm_TxN(1,i) < epsilon
            B(:,i) = zeros(3, 1);
            k(:,i) = zeros(3, 1);
        else
            B(:,i) = TxN(:,i) / norm_TxN(1,i);
            k(:,i) = TxN(:,i) / (norm_TxN(1,i)^3);
        end
        
        
        t(i) = -B(:,i)'*N(:,i);
    end
    
end

%% plot O/P results
sizeCRF = 1;
lineWidth = 1;
if metricUnitType==1  
    labelMetricUnit= 'm';
elseif metricUnitType==0
    labelMetricUnit= 'mm';
    sizeCRF=10;
    lineWidth = 3;
end
set(0,'defaulttextinterpreter','latex')

%figure(1), 
figure('Color',[1 1 1]), clf, hold on;
plot3(w_refPathCurve(1,:), w_refPathCurve(2,:), w_refPathCurve(3,:), 'k--.')%, 'LineWidth', 1.5);
view(3),
grid on, 
%hold off;

xlabel(['$x_{w}$(',labelMetricUnit,')'], 'fontweight','bold');  
ylabel(['$y_{w}$(',labelMetricUnit,')'], 'fontweight','bold');  
zlabel(['$z_{w}$(',labelMetricUnit,')'], 'fontweight','bold');
%BLegend = legend('desired path', 'approach path', 'insertion path');
%set(BLegend, 'fontsize', 14, 'location', 'NorthEast');

for i=1:3:lengthPath
    quiver3(w_refPathCurve(1,i), w_refPathCurve(2,i), w_refPathCurve(3,i), sizeCRF*T(1,i), sizeCRF*T(2,i), sizeCRF*T(3,i), 'color', 'g', 'LineWidth', lineWidth);
    quiver3(w_refPathCurve(1,i), w_refPathCurve(2,i), w_refPathCurve(3,i), sizeCRF*N(1,i), sizeCRF*N(2,i), sizeCRF*N(3,i), 'color', 'r', 'LineWidth', lineWidth);
    quiver3(w_refPathCurve(1,i), w_refPathCurve(2,i), w_refPathCurve(3,i), sizeCRF*B(1,i), sizeCRF*B(2,i), sizeCRF*B(3,i), 'color', 'b', 'LineWidth', lineWidth);
    %pause(0.25)    
end
hold off;

%% generate the animated GIF
if op_animatedGIF    
    
    %filename = 'testAnimated.gif'; % Specify the output file name
    if curveType == 0
        filename = 'spiralCurve.gif';
        writerObj_forAVI = VideoWriter('spiralCurve.avi');
        %saveas(gcf,'spiralCurveComplete.png')
        printeps(1,'spiralCurveComplete');
        opDir_forPNG = 'opImages_spiralCurve';
    elseif curveType == 1
        filename = 'zigzagCurve.gif';
        writerObj_forAVI = VideoWriter('zigzagCurve.avi');
        printeps(1,'zigzagCurveComplete');
        opDir_forPNG = 'opImages_zigzagCurve';
    elseif curveType == 2
        filename = 'urocsCurve.gif';
        writerObj_forAVI = VideoWriter('urocsCurve.avi');
        printeps(1,'urocsCurveComplete');
        opDir_forPNG = 'opImages_urocsCurve';
    elseif curveType == 3
        filename = 'curvedCurve.gif';
        writerObj_forAVI = VideoWriter('curvedCurve.avi');
        printeps(1,'curvedCurveComplete');
        opDir_forPNG = 'opImages_curvedCurve';
    elseif curveType == 4
        filename = 'lissajousCurve.gif';
        writerObj_forAVI = VideoWriter('lissajousCurve.avi');
        printeps(1,'lissajousCurveComplete');
        opDir_forPNG = 'opImages_lissajousCurve';
    end
    
    % check if the folder of PNG images exists
    %opDir_forPNG = 'opImages';
    flagExistFolder = exist(opDir_forPNG, 'dir');
    if flagExistFolder~=7 && op_saveImagesPNG
        mkdir(opDir_forPNG);
    end
    
    % activate the grapper for the video AVI
    if op_animatedAVI
        writerObj_forAVI.Quality=100;
        writerObj_forAVI.FrameRate=5;
        open(writerObj_forAVI); 
    end
    
    timeDelay = 0.1;
    figureToSave = figure('Color',[1 1 1]);
    for i=1:lengthPath
        figure(figureToSave), clf, hold on;
        plot3(w_refPathCurve(1,:), w_refPathCurve(2,:), w_refPathCurve(3,:), 'k--.');%, 'LineWidth', 1.5);

        quiver3(w_refPathCurve(1,i), w_refPathCurve(2,i), w_refPathCurve(3,i), sizeCRF*T(1,i), sizeCRF*T(2,i), sizeCRF*T(3,i), 'color', 'g', 'LineWidth', lineWidth);
        quiver3(w_refPathCurve(1,i), w_refPathCurve(2,i), w_refPathCurve(3,i), sizeCRF*N(1,i), sizeCRF*N(2,i), sizeCRF*N(3,i), 'color', 'r', 'LineWidth', lineWidth);
        quiver3(w_refPathCurve(1,i), w_refPathCurve(2,i), w_refPathCurve(3,i), sizeCRF*B(1,i), sizeCRF*B(2,i), sizeCRF*B(3,i), 'color', 'b', 'LineWidth', lineWidth);
        

        xlabel(['$x_{w}$(',labelMetricUnit,')'], 'fontweight','bold');  
        ylabel(['$y_{w}$(',labelMetricUnit,')'], 'fontweight','bold');  
        zlabel(['$z_{w}$(',labelMetricUnit,')'], 'fontweight','bold');
        %BLegend = legend('desired path', 'approach path', 'insertion path');
        %set(BLegend, 'fontsize', 14, 'location', 'NorthEast');
        
        view(3), 
        grid on, 
        axis tight, axis equal; %axis off;
        hold off;
        %drawnow

        % Capture the plot as an image
        frame = getframe(figureToSave);
        im{i} = frame2im(frame);    
        %[A,map] = rgb2ind(im{i},256);
        [A,map] = rgb2ind(frame.cdata,256,'nodither');
        
        if op_saveImagesPNG
            imwrite(A,map,[opDir_forPNG,'/im_',num2str(i),'.png']);
            %printeps(1,[opDir,'/im_',num2str(i)]);
        end

        % Write to the GIF File
        if i == 1
            imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',timeDelay);
        else
            imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',timeDelay);
        end
        
        if op_animatedAVI
            writeVideo(writerObj_forAVI,frame); 
        end
    end %for i=1:lengthPath
    
    if op_animatedAVI
        close(writerObj_forAVI) 
    end
end %if op_animatedGIF


