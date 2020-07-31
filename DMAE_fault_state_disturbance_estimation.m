%% ----------------------------------------------------------------------
% This is the code for the following journal paper:
%
%      P. Lu, E. van Kampen, C. C. de Visser, Q. P. Chu
%      Framework for state and unknown input estimation of linear time-varying systems
%      Automatica, 2016(73), 145-154
%
%   This program shows the results of Example 2, Section 5 of the paper.
%
%   Please read the paper for more details.
%
%   If you have more questions, please contact the author:
%
%   Author: Peng Lu, Delft University of Technology
%
%   email:  P.Lu-1@tudelft.nl
%
%   released in June 2016
%----------------------------------------------------------------------

%% ---------------- start of the program --------------------------

close all
clear all


% choose three cases, refer to the paper for more details
disp('0. only disturbances, without faults'); %
disp('1. no disturbances, only faults'); %
disp('2. both disturbances and faults'); %

f_flag = input('Choose a number (0,1,2): ');


% load the corresponding system data
% these data are generated using "generate_sys_data.m"
if f_flag == 0
    load sys_disturbance.mat
elseif f_flag == 1
    load sys_fault.mat
elseif f_flag == 2
    load sys_disturbance_fault.mat
else
    disp('I told you to select from 0,1,2');
end


% sampling rate
delta_t = 1;

% initial
x_ob_0  =  zeros(2,1);
P_x0  =  1e0*eye(2);


% initial condition for the two filters
mean_x0_ff=[x_ob_0; 0; 0; ];

mean_x0_M1=[x_ob_0; 0; 0; 0; 0;];

%
[dim_sys,~]   =   size(mean_x0_ff);
[dim_out,~]   =   size(z_real);
[dim_dis,~]  =  size(dk);
[dim_f,~]  =  size(fk);

P_x0_ff=1e0*eye(dim_sys);
P_x0_M1=1e0*eye(dim_sys+dim_f);

% initialization of matrices
x_ob   =   zeros(dim_sys,2*dim_sys+1);
z_ob   =   zeros(dim_out,2*dim_sys+1);
x_ob_mean   =   zeros(dim_sys,gen);
z_ob_mean   =   zeros(dim_out,gen);
inno   =   zeros(dim_out,gen);
norm_inno   =   zeros(1,gen);
error_x   =   zeros(2,gen);
norm_error_x   =   zeros(dim_sys,gen);
x_ob_filter   =   zeros(dim_sys,gen);
z_ob_filter   =   zeros(dim_out,gen);
P_ob_filter   =   zeros(dim_sys,dim_sys,gen);
residual   =   zeros(dim_out,gen);
norm_residual   =   zeros(1,gen);

d_ob   =   zeros(dim_dis,gen);
Pd_ob   =   zeros(dim_dis,dim_dis,gen);

% noise covariance
Qk  =  Qk_real;
Rk  =  Rk_real;



%--------------  reconstruct the system matrices ------------
%-------- system matrices for the fault-free filter
Phi = [                Phi                 E;
        zeros(dim_dis,dim_sys-dim_dis) eye(dim_dis) ];
  
B = [B; zeros(dim_dis,1)];
H = [H zeros(dim_out,dim_dis)];    
%--------------  

% initial for MM
[dim_sys_MM,~]=size(mean_x0_M1);

inno_ff=zeros(dim_out,gen);
error_x_ff=zeros(2,gen);
x_ob_filter_ff=zeros(dim_sys,gen);
Pzk_ff=zeros(dim_out,dim_out,gen);
P_ob_filter_ff=zeros(dim_sys,dim_sys,gen);

inno_M1=zeros(dim_out,gen);
error_x_M1=zeros(2,gen);
x_ob_filter_M1=zeros(dim_sys_MM,gen);
Pzk_M1=zeros(dim_out,dim_out,gen);
P_ob_filter_M1=zeros(dim_sys_MM,dim_sys_MM,gen);

%
x_ob_filter_pdf = zeros(dim_sys,gen);
x_ob_filter_pdf2 = zeros(dim_sys,gen);
d_ob_ff = zeros(dim_dis,gen);
d_ob_M1 = zeros(dim_dis,gen);

pd_ff=zeros(1,gen);
pd_M1=zeros(1,gen);


p_t=zeros(1,gen);
p_ff=zeros(1,gen);
p_M1=zeros(1,gen);

p_k=zeros(2,gen);

liho_ff=zeros(1,gen);
liho_M1=zeros(1,gen);


d_flag_gen = zeros(1,gen);

S_inno_ff=zeros(dim_out,dim_out,gen);
C_inno_ff=zeros(dim_out,dim_out,gen);
S_inno_M1=zeros(dim_out,dim_out,gen);
C_inno_M1=zeros(dim_out,dim_out,gen);
norm_S_inno_ff = zeros(1,gen); 
norm_S_inno_M1 = zeros(1,gen); 
norm_C_inno_ff = zeros(1,gen); 
norm_C_inno_M1 = zeros(1,gen); 
Fading_MM = zeros(dim_out,dim_out,gen);
fading_gen = zeros(1,gen); 
Qd_est_ff = zeros(dim_sys,dim_sys,gen);
Qd_est_MM = zeros(dim_sys_MM,dim_sys_MM,gen);

Fading_ff = eye(4);
Fading_M1 = eye(6);

K_ff=zeros(dim_sys,dim_out,gen);
K_M1=zeros(dim_sys_MM,dim_out,gen);


% pdf min & max
I_min=zeros(1,gen);
I_max=zeros(1,gen);

% the initial PDF
p0_0=0.95;
p1_0=0.05;

g = 9.81;
 
% true model probabilities
PI=zeros(2,gen);
if f_flag ~= 0
    for k=1:500
        if norm(fk(:,k)) > 1e-6
            PI(2,k) = 1;
            PI(1,k) = 0;
        else
            PI(2,k) = 0;
            PI(1,k) = 1;
        end
    end
else
    PI(1,1:500) = 1;
    PI(2,1:500) = 0;
end




% parameters 
d_flag = 0;


% the inverse of H
inv_H = eye(2)/H(1:2,1:2);

i_max = 1; % initial value

%---------------------------- run the DMAE ----------------------------
for k=1:gen
    
    % augmentation
    [inno_M1(:,k),error_x_M1(:,k),K_M1(:,:,k),Pzk_M1(:,:,k),x_ob_filter_M1(:,k),P_ob_filter_M1(:,:,k)]...
        = EKF_fault_filter(f_flag,i_max,Phi,B,H,F,mean_x0_M1,P_x0_M1,x_real,z_real,Qk,Rk,u,dim_sys_MM,dim_dis,dim_f,k,Qd_est_MM(:,:,k));

    % no augmentation
    [inno_ff(:,k),error_x_ff(:,k),K_ff(:,:,k),Pzk_ff(:,:,k),x_ob_filter_ff(:,k),P_ob_filter_ff(:,:,k)]...
        = EKF_fault_free_filter(Phi,B,H,mean_x0_ff,P_x0_ff,x_real,z_real,Qk,Rk,u,dim_sys,dim_dis,k,Qd_est_ff(:,:,k));

    %------------------------------ PDF update -----------------------
    % likelihood function
    liho_ff(:,k)=(inno_ff(:,k)'/(2*Pzk_ff(:,:,k))*inno_ff(:,k));
    liho_M1(:,k)=(inno_M1(:,k)'/(2*Pzk_M1(:,:,k))*inno_M1(:,k));

    pd_ff(:,k)=exp(-liho_ff(:,k))/(det(2*pi*Pzk_ff(:,:,k))^(1/2));
    pd_M1(:,k)=exp(-liho_M1(:,k))/(det(2*pi*Pzk_M1(:,:,k))^(1/2));

    pd0=pd_ff(:,k)*p0_0;
    p1=pd_M1(:,k)*p1_0;

    % sum of the PDF
    p_t(:,k)=pd0+p1;
    p_ff(:,k)=pd0/p_t(:,k);
    p_M1(:,k)=p1/p_t(:,k);

    
    %--------------       next gen     ---------------
    mean_x0_ff=x_ob_filter_ff(:,k);
    P_x0_ff=P_ob_filter_ff(:,:,k);

    mean_x0_M1=x_ob_filter_M1(:,k);
    P_x0_M1=P_ob_filter_M1(:,:,k);
    
    
    % mark the index for max & min pdf
    p_k(:,k)=[p_ff(:,k);p_M1(:,k)];
    [min_p,I_min(k)]=min(p_k(:,k));
    [max_p,I_max(k)]=max(p_k(:,k));
 
    
    flag_f=0;  

    %-------------- start of Selective Reinitialization algorithm --------------
    % read the paper and related papers for more details
    
    % ------------  no faults:   reinitialise fault filter using fault-free filter  ------------
    if I_max(k)==1  %
 
        mean_x0_M1(1:dim_sys,:)=mean_x0_ff;  P_x0_M1(1:dim_sys,1:dim_sys)=P_x0_ff;
        mean_x0_M1(dim_sys+1:dim_sys_MM)=1e-3*ones(dim_f,1);  
        
        P_x0_M1(1:dim_sys,dim_sys+1:dim_sys_MM)=zeros(dim_sys,dim_f);
        P_x0_M1(dim_sys+1:dim_sys_MM,1:dim_sys)=zeros(dim_f,dim_sys);       
        P_x0_M1(dim_sys+1:dim_sys_MM,dim_sys+1:dim_sys_MM)=1e2*eye(dim_f);
        i_max = I_max(k);
    end
    % -----------------------------------    end    -------------------------------

    % --------------  faults:   reinitialise fault-free filter using fault filter    ------------
    if k>1 && I_max(k)==2  

        mean_x0_ff(1:dim_sys,:)=mean_x0_M1(1:dim_sys,:);
        P_x0_ff(1:dim_sys,1:dim_sys)=P_x0_M1(1:dim_sys,1:dim_sys);
        i_max = I_max(k);
    end
    % -------------------------------    end    -------------------------------------
   
    %-------------- end of Selective Reinitialization algorithm --------------
    
    % min, prevent lock out
    p_ff(:,k)=max(p_ff(:,k),0.001);
    p_M1(:,k)=max(p_M1(:,k),0.001);

    % max, prevent lock out
    p_ff(:,k)=min(p_ff(:,k),0.999);
    p_M1(:,k)=min(p_M1(:,k),0.999);
    
    % the next 
    p0_0=p_ff(:,k);
    p1_0=p_M1(:,k);
    
    
 
    %----------------------- fusion of estimates ------------------------
    x_ob_filter_pdf(:,k)=p_ff(:,k)*x_ob_filter_ff(:,k)+p_M1(:,k)*x_ob_filter_M1(1:dim_sys,k);
    
    % choose the model with the highest PDF, do not fuse all the models
    [p_max,I]=max([p_ff(:,k),p_M1(:,k)]);
    x_com=[x_ob_filter_ff(:,k),x_ob_filter_M1(1:dim_sys,k)];
    x_ob_filter_pdf2(:,k)=x_com(:,I);

    
    %--------------- adaptive estimation of Qd ------------------
    % read the paper and related papers for more details
    %----------------        innovation     ----------------------
    S_inno_ff(:,:,k) = inno_ff(:,k)*inno_ff(:,k)';
    S_inno_M1(:,:,k) = inno_M1(:,k)*inno_M1(:,k)';

    %----- remove the increase caused by the reinitialization ------
    if S_inno_M1(1,1,k)>5 && k>100
        S_inno_M1(1,1,k) = S_inno_M1(1,1,k-1);
    end
    if S_inno_M1(2,2,k)>5 && k>100
        S_inno_M1(2,2,k) = S_inno_M1(2,2,k-1);
    end
    %----- remove the increase caused by the reinitialization ------
    
    N_win = 10;% width of moving window
    % sum of the innovation square
    if N_win == 1
        C_inno_ff(:,:,k) = S_inno_ff(:,:,k);
        C_inno_M1(:,:,k) = S_inno_M1(:,:,k);
    end
    if N_win > 1
        if k == 1
            C_inno_ff(:,:,k) = S_inno_ff(:,:,k);
            C_inno_M1(:,:,k) = S_inno_M1(:,:,k);
        elseif k >1 && k < N_win
            C_inno_ff(:,:,k) = sum(S_inno_ff(:,:,1:k),3)/(k-1);
            C_inno_M1(:,:,k) = sum(S_inno_M1(:,:,1:k),3)/(k-1);
        elseif k > N_win
            C_inno_ff(:,:,k) = sum(S_inno_ff(:,:,k-N_win+1:k),3)/(N_win-1);
            C_inno_M1(:,:,k) = sum(S_inno_M1(:,:,k-N_win+1:k),3)/(N_win-1);
        end 
    end
    norm_S_inno_ff(:,k) = norm(S_inno_ff(:,:,k));
    norm_S_inno_M1(:,k) = norm(S_inno_M1(:,:,k));
    norm_C_inno_ff(:,k) = norm(C_inno_ff(:,:,k));
    norm_C_inno_M1(:,k) = norm(C_inno_M1(:,:,k));
    %----------------        innovation     ----------------------
 
    % use Fault Model to est
    Qd_0 = C_inno_M1(:,:,k) - H(1:2,1:2)*Qk*H(1:2,1:2)' - Rk;
    
    % can guarantee positive
    Qd_est = diag([ max(Qd_0(1,1),0), max(Qd_0(2,2),0)]);
    %
    Qdx_est = inv_H*Qd_est*inv_H';
    %--------------------- take account of the E matrix ---------------------

    if E == zeros(2,2) 
        Qdd_est = Qdx_est;
    else
        Qdd_est = inv(E) *Qdx_est* inv(E'); % Qd
    end
    Qd_est_ff(:,:,k+1) = diag([ Qdx_est(1,1), Qdx_est(2,2), Qdd_est(1,1), Qdd_est(2,2) ]);
    Qd_est_MM(:,:,k+1) = diag([ Qdx_est(1,1), Qdx_est(2,2), Qdd_est(1,1), Qdd_est(2,2), 0, 0 ]);
    
    %----------- take the correlation into account --------------
    % this is optional
    %----------- see the paper why we add this
    E_Qd = E * Qdd_est;
    Qd_E = Qdd_est * E';
    Qd_est_ff(:,:,k+1) = zeros(4,4);
    Qd_est_MM(:,:,k+1) = zeros(6,6);
    
    Qd_est_ff(1:2,1:2,k+1) = Qdx_est + Qk;
    Qd_est_ff(1:2,3:4,k+1) = E_Qd; % can also be set to zero
    Qd_est_ff(3:4,1:2,k+1) = Qd_E; % can also be set to zero
    Qd_est_ff(3:4,3:4,k+1) = Qdd_est;
    Qd_est_MM(1:4,1:4,k+1) = Qd_est_ff(1:4,1:4,k+1);




    
end


%
Time=delta_t*(1:k);


% RMSE root mean square error 
RMSE=sqrt(mean(error_x_ff.^2,2));

RMSE_M1=sqrt(mean(error_x_M1.^2,2));
fprintf(1,'the RMSE is:\n%f\n%f\n%f\n ',sqrt(mean((error_x_M1(1,1:k)).^2)),sqrt(mean((error_x_M1(2,1:k)).^2)));
fprintf(1,'the norm of root mean square error is:\n%f\n\n ',norm(RMSE_M1(:)));



% plots
%--------------- model probabilities
figure;hold on;axis off;
subplot(211);hold on; plot(Time,PI(1,1:k),'r','linewidth',2); plot(Time,p_ff(1:k),'b--','linewidth',2); ylabel('p_{nf}','fontsize',12);set(gca,'ylim',[0 1],'fontsize',12);
% box off; set(gca,'position',[0.1 0.76 0.85 0.15],'fontsize',12);
subplot(212);hold on; plot(Time,PI(2,1:k),'r','linewidth',2); plot(Time,p_M1(1:k),'b--','linewidth',2); ylabel('p_{af}','fontsize',12);set(gca,'ylim',[0 1],'fontsize',12);
h=legend('True','DMAE'); set(h,'color','none','edgecolor','white');
h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time (s)','fontsize',12)
set(h1,'Box','off')



%--------------- state estimation 
figure;
subplot(211); hold on; plot(x_real(1,1:k),'r','linewidth',2); plot(x_ob_filter_pdf2(1,1:k),'b--','linewidth',2); ylabel('x_1','fontsize',12);grid; set(gca,'fontsize',12);
subplot(212); hold on; plot(x_real(2,1:k),'r','linewidth',2); plot(x_ob_filter_pdf2(2,1:k),'b--','linewidth',2); ylabel('x_2','fontsize',12);grid; set(gca,'fontsize',12);
h=legend('True','DMAE'); set(h,'color','none','edgecolor','white');
h2=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time (s)','fontsize',12)
set(h2,'Box','off')



%--------------- disturbance estimation 
figure;
subplot(211); hold on; plot(Time,dk(1,1:k),'r','linewidth',2); plot(Time,x_ob_filter_pdf2(3,1:k),'b--','linewidth',2);  ylabel('d_1 ','fontsize',12);grid; set(gca,'fontsize',12);% set(gca,'xlim',[0 delta_t*k],'ylim',[-0.2 1.2],'fontsize',12);
subplot(212); hold on; plot(Time,dk(2,1:k),'r','linewidth',2); plot(Time,x_ob_filter_pdf2(4,1:k),'b--','linewidth',2); ylabel('d_2','fontsize',12);grid; set(gca,'fontsize',12);%  set(gca,'xlim',[0 delta_t*k],'ylim',[-0.02 0.1],'fontsize',12);
h2=legend('True','DMAE'); set(h2,'color','none','edgecolor','white');
h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time (s)','fontsize',12)
set(h1,'Box','off')

error_d1 = dk(1,1:k) - x_ob_filter_pdf2(3,1:k);
error_d2 = dk(2,1:k) - x_ob_filter_pdf2(4,1:k);
fprintf(1,'the RMSE of the disturbance estimation is:\n %f\n %f\n ',sqrt(mean((error_d1).^2)),sqrt(mean((error_d2).^2)));


%--------------- fault estimation 
figure;
subplot(211); hold on; plot(Time,fk(1,1:k),'r','linewidth',2); plot(Time,p_M1(1:k).*x_ob_filter_M1(5,1:k),'b--','linewidth',2);  ylabel('f_1','fontsize',12);grid; set(gca,'fontsize',12);% set(gca,'xlim',[0 delta_t*k],'ylim',[-0.2 1.2],'fontsize',12);
subplot(212); hold on; plot(Time,fk(2,1:k),'r','linewidth',2); plot(Time,p_M1(1:k).*x_ob_filter_M1(6,1:k),'b--','linewidth',2); ylabel('f_2','fontsize',12);grid; set(gca,'fontsize',12);%  set(gca,'xlim',[0 delta_t*k],'ylim',[-0.02 0.1],'fontsize',12);
h2=legend('True','DMAE'); set(h2,'color','none','edgecolor','white');
h1=axes('position',[0.52 0.0001 0.0001 0.0001],'fontsize',12); title('time (s)','fontsize',12)


%
error_f1 = fk(1,1:k) - p_M1(1:k).*x_ob_filter_M1(5,1:k);
error_f2 = fk(2,1:k) - p_M1(1:k).*x_ob_filter_M1(6,1:k);
fprintf(1,'the RMSE of the fault estimation is:\n %f\n %f\n ',sqrt(mean((error_f1).^2)),sqrt(mean((error_f2).^2)));




