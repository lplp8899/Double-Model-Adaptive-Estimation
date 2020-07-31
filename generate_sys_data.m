%% ----------------------------------------------------------------------
% This is the code for the following journal paper:
%
%      P. Lu, E. van Kampen, C. C. de Visser, Q. P. Chu
%      Framework for state and unknown input estimation of linear time-varying systems
%      Automatica, 2016
%
%   This program generates the data used in Example 2, Section 5 of the paper.
%
%   If you have more questions, please contact the author:
%
%   Author: Peng Lu, Delft University of Technology
%
%   email:  P.Lu-1@tudelft.nl
%
%   released in June 2016
%----------------------------------------------------------------------

%% -------------------- start of the program --------------------------
% this program generates the data using the system model
close all
clear all
% clc

% choose three cases, refer to the paper for more details
disp('0. only disturbances, without faults'); %
disp('1. no disturbances, only faults'); %
disp('2. both disturbances and faults'); %

f_flag = input('Choose: ');

% load the disturbance data
load dist.mat


%
gen = 500; % time steps
dim_sys = 2;
dim_out = 2;
dim_d = 2;
x_real=zeros(dim_sys,gen);

z_real=zeros(dim_out,gen);

%----------- system matrices -------------
% Phi
Phi = [-0.0005 -0.0084;
        0.0517  0.8069]; % it is Matrix A in the paper

B = [ 0.1815;
      1.7902];   

H = [1  0;
    0   1;];  

%
if f_flag == 0
    E = [0.629 0;
         0     -0.52504;];
    F = [0  0;
         0  0;];
elseif f_flag == 1
    E = [0 0;
         0 0;];
    F = [1  0;
         0  1;]; 
elseif f_flag == 2
    E = [0.629 0;
        0     -0.52504;];
    F = [1  0;
         0  1;];
end

Qk_real = diag([0.002^2,0.002^2]);
Rk_real = diag([0.01^2,0.01^2]);

dt = 1;

% generate random noise
randn('state',1000)
w = sqrt(Qk_real)*randn(2,gen);
randn('state',2000)
v = sqrt(Rk_real)*randn(2,gen);

% dk
dk = zeros(dim_d,gen);
if f_flag ~= 1
    dk(1,1:500) = d1(1,1:500);
    dk(2,1:500) = d2(1,1:500);
end

% fk
fk = zeros(dim_d,gen);
if f_flag ~= 0
    fk(1,100:200) = 4;
    fk(2,300:400) = 3;
end

% known input
u = zeros(1,gen);
u(1:200) = 0.5;
u(201:300) = -0.5;
u(301:500) = 0.5;

% initial state
x_real_0 = zeros(2,1);

%---------------------- generate system data -----------------------
for k=1:gen
    x_real(:,k) = Phi*x_real_0 + B*u(:,k) + E*dk(:,k) + w(:,k);
    z_real(:,k) = H*x_real(:,k) + F*fk(:,k) + v(:,k);
    x_real_0 = x_real(:,k);
end
Time = dt*(1:k);



% plots
figure;
subplot(211); hold on; grid; plot(Time,x_real(1,1:k),'b','linewidth',1.5); 
subplot(212); hold on; grid; plot(Time,x_real(2,1:k),'b','linewidth',1.5); 


figure;
subplot(211); hold on; grid; plot(Time,z_real(1,1:k),'b','linewidth',1.5); plot(Time,x_real(1,1:k),'r','linewidth',1.5);
subplot(212); hold on; grid; plot(Time,z_real(2,1:k),'b','linewidth',1.5); plot(Time,x_real(2,1:k),'r','linewidth',1.5);
legend('measurement','state','fontsize',13);


figure;
subplot(211); hold on; plot(Time,dk(1,1:k),'r','linewidth',2); ylabel('d1','fontsize',12);grid; set(gca,'fontsize',12);% set(gca,'xlim',[0 delta_t*k],'ylim',[-0.2 1.2],'fontsize',12);
subplot(212); hold on; plot(Time,dk(2,1:k),'r','linewidth',2); ylabel('d2','fontsize',12);grid; set(gca,'fontsize',12);%  set(gca,'xlim',[0 delta_t*k],'ylim',[-0.02 0.1],'fontsize',12);
h2=legend('True','DMAE1'); set(h2,'color','none','edgecolor','white');
h1=axes('position',[0.10 0.05 0.84 0.86],'fontsize',12);axis off;title(h1,'Disturbance')
h1=axes('position',[0.50 0.0001 0.0001 0.0001],'fontsize',12); title('time (s)','fontsize',12)


figure;
subplot(211); hold on; plot(Time,fk(1,1:k),'r','linewidth',2); ylabel('f1','fontsize',12);grid; set(gca,'fontsize',12);% set(gca,'xlim',[0 delta_t*k],'ylim',[-0.2 1.2],'fontsize',12);
subplot(212); hold on; plot(Time,fk(2,1:k),'r','linewidth',2); ylabel('f2','fontsize',12);grid; set(gca,'fontsize',12);%  set(gca,'xlim',[0 delta_t*k],'ylim',[-0.02 0.1],'fontsize',12);
h2=legend('True','DMAE1'); set(h2,'color','none','edgecolor','white');
h1=axes('position',[0.10 0.05 0.84 0.86],'fontsize',12);axis off;title(h1,'Fault ')
h1=axes('position',[0.50 0.0001 0.0001 0.0001],'fontsize',12); title('time (s)','fontsize',12)




% save the data and use it for fault and disturbance estimation
if f_flag == 0
    save sys_disturbance.mat
elseif f_flag == 1
    save sys_fault.mat
elseif f_flag == 2
    save sys_disturbance_fault.mat
end



