
clc; clear all; close all; warning off; pause(0.01)

load IODATA.mat
load x_Initial.mat

%%%%% OKID 過程 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[G_ok,H_ok,C_ok,D_ok,L_ok] = Auxi_OKID_JXL(u,y,2,3,0);

%%%%% 重現過程 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[th,xh,yh,eh,Gd,Hd,Ld,x_Initial] = Auxi_OKID_Process_step2(u,y,G_ok,H_ok,C_ok,D_ok,L_ok,(10/101),0,1e8,x_Initial);

figure(1); plot(1:length,y(1,:),'.b-',1:length,yh(1,:),'xr:'); xlim([0,length]); 
ylabel('y_1 vs. yh_1');
xlabel('Time (sec)');

figure(2); plot(1:length,y(1,:)-yh(1,:),'k'); xlim([0,length]); 
ylabel('y_1 － yh_1');
xlabel('Time (sec)');

figure(3); plot(1:length,y(2,:),'.b-',1:length,yh(2,:),'xr:'); xlim([0,length]); 
ylabel('y_2 vs. yh_2');
xlabel('Time (sec)');

figure(4); plot(1:length,y(2,:)-yh(2,:),'k'); xlim([0,length]); 
ylabel('y_2 － yh_2');
xlabel('Time (sec)');

save x_Initial.mat x_Initial