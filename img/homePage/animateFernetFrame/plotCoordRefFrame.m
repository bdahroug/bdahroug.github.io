function [  ] = plotCoordRefFrame( i_M_j, originName, sizeCRF)

    set(0,'defaulttextinterpreter','latex')
    if (nargin<3) 
        sizeCRF=1;
    end
    if (nargin<2) 
        originName='o';
        sizeCRF=1;
    end
    
    hold_state=true;
    if (~ishold)
        view(3), hold on;
        hold_state=false;
    end
    
    %the position of origion
    i_P_j= i_M_j(1:3,4);
    
    u0=[1;0;0;0];  v0=[0;1;0;0];  w0=[0;0;1;0];
    u1=i_M_j*u0;   v1=i_M_j*v0;   w1=i_M_j*w0;
    u1=sizeCRF*u1;   v1=sizeCRF*v1;   w1=sizeCRF*w1;
    quiver3(i_P_j(1),i_P_j(2),i_P_j(3),u1(1),u1(2),u1(3),'r','linewidth',3);
    quiver3(i_P_j(1),i_P_j(2),i_P_j(3),v1(1),v1(2),v1(3),'g','linewidth',3);
    quiver3(i_P_j(1),i_P_j(2),i_P_j(3),w1(1),w1(2),w1(3),'b','linewidth',3);
    
    fontSize=10;
    text(i_P_j(1),i_P_j(2),i_P_j(3), originName,'Color','k','Fontsize',fontSize*2);
    text(i_P_j(1)+u1(1),i_P_j(2)+u1(2),i_P_j(3)+u1(3), ['$x_{',originName,'}$'],'interpreter','latex', 'Color', 'r', 'Fontsize',fontSize,'Fontweight','bold');
    text(i_P_j(1)+v1(1),i_P_j(2)+v1(2),i_P_j(3)+v1(3), ['$y_{',originName,'}$'],'interpreter','latex', 'Color', 'g', 'Fontsize',fontSize,'Fontweight','bold');
    text(i_P_j(1)+w1(1),i_P_j(2)+w1(2),i_P_j(3)+w1(3), ['$z_{',originName,'}$'],'interpreter','latex', 'Color', 'b', 'Fontsize',fontSize,'Fontweight','bold');


    if (~hold_state) 
        hold off
    end
    
end

