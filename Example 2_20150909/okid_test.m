clc; clear all; close all;
%% System matrix
A =  [-1.0, 0.0, 0.00,-0.5;     % ╰参把计A
      0.0,-0.5,-0.25,-0.5;
      0.0, 0.0,-0.50, 0.0;
      0.0, 0.0, 0.00,-0.5]; 

B =[ 0.50, 0.50 1.0;            % ╰参把计B
     -0.25,-0.25 0.7;
      0.50, 0.50 -0.9;
      0.50,-0.50 0.5];

 
C = [ 1, 1, 0,-1.5;             % ╰参把计C
      0, 1, 0,-1.0];
 
D = [0, 0 ,0;                   % ╰参把计D
     0, 0, 0];
 
 initial_state=[0.1             % 篈﹍
                0.2
                0.3
                0.4];

[p m] = size(D);                % p块蝴m块蝴
n = size(A);                    % n篈蝴

Ts = 0.1;                       % 瞒床蹦妓丁
Tend = 10;
[G,H] = c2d(A,B,Ts);            % 盢╰参把计暗瞒床て

%% Reference
t_ds = 0:Ts:Tend;
u = 0.2*randn(m,size(t_ds,2));

ii = 0;
x(:,1) = initial_state;
for t = 1 : size(t_ds,2)
    ii = ii+1;
    x(:,ii+1) = G*x(:,ii) + H*u(:,ii);
    y(:,ii) = C*x(:,ii);
end
length = size(u,2);
save IODATA.mat u y length