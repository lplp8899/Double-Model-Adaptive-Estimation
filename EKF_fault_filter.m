%----------------------------------------------------------------------
%   Author: Peng Lu, Delft University of Technology 
%   email:  P.Lu-1@tudelft.nl
%   released in June 2016
%----------------------------------------------------------------------

function [inno,error_x,K,V,x_ob_filter,P_ob_filter]...
    = EKF_fault_filter(f_flag,i_max,Phi,B,H,F,x_ob_0,P_x0,x_real,z_real,Qk,Rk,u,dim_sys,dim_dis,dim_f,k,Qd_est_MM)


if f_flag == 2 && i_max == 2
    
    % refer to section 4.3 of the paper
    x_ob(1:4,:) = Phi*x_ob_0(1:4) + B*u(:,k);

    Ga = [eye(2); zeros(dim_dis,2)];
    P = Phi*P_x0(1:4,1:4)*Phi' + Ga*Qk*Ga';
    Qkk = 0;


    P = P + Qkk + Qd_est_MM(1:4,1:4);

    z_ob = H*x_ob;
    
    % residual
    inno =  z_real(:,k)-z_ob - x_ob_0(5:6,:);
 
    V=H*P*H'+Rk;

    K(1:4,:)=P*H'/V;
    
    K(5:6,:) = zeros(2,2);

    %predict    
    x_ob_filter(1:4,:) = x_ob+K(1:4,1:2)*inno;
    P_ob_filter(1:4,1:4,:) =(eye(dim_sys-2)-K(1:4,1:2)*H)*P*(eye(dim_sys-2)-K(1:4,1:2)*H)'...
        +K(1:4,1:2)*Rk*K(1:4,1:2)';

    %
    x_ob_filter(5:6,:) = x_ob_0(5:6,:); 
    P_ob_filter(1:4,5:6) = zeros(4,2);
    P_ob_filter(5:6,1:4) = zeros(2,4);
    P_ob_filter(5:6,5:6) = P_x0(5:6,5:6);
    
    %error & residual
    error_x=x_real(:,k)-x_ob_filter(1:2,1);
    

    
else
    
%------------------      states       ------------------
%-------- system model for the fault filter
%-------- augmented based on fault-free filter
    Phi_a = [ Phi zeros(dim_sys-dim_f,dim_f);
            zeros(dim_f,dim_sys-dim_f) eye(dim_f)];
    B_a = [B;zeros(dim_f,1);];
    Ga_a = [eye(2); zeros(dim_f+dim_dis,2)]; % because E_d is augmented


    x_ob = Phi_a*x_ob_0 + B_a*u(:,k) ;

    P = Phi_a*P_x0*Phi_a' + Ga_a*Qk*Ga_a';

    P = P + Qd_est_MM;


    H_a = [H F];

    z_ob = H_a*x_ob;
    
    % residual
    inno=z_real(:,k)-z_ob;


    V=H_a*P*H_a'+Rk;

    K=P*H_a'/V;


    %predict    
    x_ob_filter=x_ob+K*inno;
    P_ob_filter=(eye(dim_sys)-K*H_a)*P*(eye(dim_sys)-K*H_a)'+K*Rk*K';

       
    %error & residual
    error_x=x_real(:,k)-x_ob_filter(1:2,1);






 

end
    
    
    