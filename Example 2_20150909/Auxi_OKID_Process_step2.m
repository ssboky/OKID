% Ver. �G20130926
% Desc.�G�N OKID �ұo�� G,H,C,D,Lo�A�åH�u�D�w�����v�Ρu�w�����v���覡�i�歡�N��A�D�o�� x,y
% Info.�Gu(��ڿ�J),y(��ڿ�X), G,H,C,D,Lo, Ts(�ļˮɶ�), Predict_OKID(���{�覡)
function [th,xh,yh,e,Gd,Hd,Lo,x_in] = Auxi_OKID_Process_step2(u,y,G,H,C,D,Lo,Ts,Predict_OKID,Qo,x_Initial)
% �d�d �����Ѽƹw�]�ȳ]�w �d�d
Check_POKID = exist('Predict_OKID');
if (Check_POKID == 0 | Predict_OKID == 0)
    Predict_OKID = 0;
end
Check_Qo = exist('Qo');
if (Check_Qo == 0)
    Qo = 1e6;
end
[n,m] = size(H);
[p,n] = size(C);
% �d�d �w�����t�ΡB��J�M�[�����ѼƳ]�w �d�d    
if Predict_OKID == 1
    [A,B] = d2c(G,H,Ts);
    Qo = Qo*eye(n);
    Ro = eye(p);
    [Lc_temp,Po] = lqr(A',C',Qo,Ro);
    Lc = Lc_temp';
    Lo = (G-eye(n))*inv(A)*Lc*inv(eye(p)+C*(G-eye(n))*inv(A)*Lc);   
end
Gd = (G-Lo*C*G);
Hd = (H-Lo*C*H);
% �d�d ��l�ȳ]�w �d�d
Num_Sample = length(u);
xh(:,1) = x_Initial;
yh(:,1) = C*xh(:,1)+D*u(:,1);
e(:,1) = y(:,1)-yh(:,1);
% �d�d �����L�{ �d�d
for i = 2:Num_Sample
    if Predict_OKID == 0 % �D�w����
        xh(:,i) = G*xh(:,i-1)+H*u(:,i-1)-Lo*e(:,i-1);
    elseif Predict_OKID == 1 % �w����
        xh(:,i) = Gd*xh(:,i-1)+Hd*u(:,i-1)+Lo*y(:,i);
    end
    yh(:,i) = C*xh(:,i)+D*u(:,i);
    e(:,i) = y(:,i)-yh(:,i);
    th(:,i) = (i-1)*Ts;
end
% �d�d ��窱�A���פε��� �d�d
x_in = xh(:,10);
disp('   Row Col')
disp(['xh',char(9),num2str(size(xh)),char(10),...
      'ue',char(9),num2str(size(u)),char(10),...
      'yh',char(9),num2str(size(yh)),char(10),...
      'th',char(9),num2str(size(th)),char(10)])