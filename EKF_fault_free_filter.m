%----------------------------------------------------------------------
%   Author: Peng Lu, Delft University of Technology 
%   email:  P.Lu-1@tudelft.nl
%   released in June 2016
%----------------------------------------------------------------------

function [inno,error_x,K,V,x_ob_filter,P_ob_filter]...
    = EKF_fault_free_filter(Phi,B,H,x_ob_0,P_x0,x_real,z_real,Qk,Rk,u,dim_sys,dim_dis,k,Qd_est_ff)

 % predict
    x_ob = Phi*x_ob_0 + B*u(:,k) ;

    Ga = [eye(2); zeros(dim_dis,2)];
    P = Phi*P_x0*Phi' + Ga*Qk*Ga';

    
    P = P + Qd_est_ff;


    z_ob = H*x_ob;
    
    % residual
    inno=z_real(:,k)-z_ob;

   
    V=H*P*H'+Rk;

    K=P*H'/V;


    %predict    
    x_ob_filter=x_ob+K*inno;
    P_ob_filter=(eye(dim_sys)-K*H)*P*(eye(dim_sys)-K*H)'+K*Rk*K';

       
    %error 
    error_x=x_real(:,k)-x_ob_filter(1:2,1);





