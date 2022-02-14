% Ver. ：20130926
% Desc.：將 OKID 所得的 G,H,C,D,Lo，並以「非預測型」或「預測型」的方式進行迭代後，求得其 x,y
% Info.：u(實際輸入),y(實際輸出), G,H,C,D,Lo, Ts(採樣時間), Predict_OKID(重現方式)
function [th,xh,yh,e,Gd,Hd,Lo,x_in] = Auxi_OKID_Process_step2(u,y,G,H,C,D,Lo,Ts,Predict_OKID,Qo,x_Initial)
% ▃▃ 相關參數預設值設定 ▃▃
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
% ▃▃ 預測型系統、輸入和觀測器參數設定 ▃▃    
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
% ▃▃ 初始值設定 ▃▃
Num_Sample = length(u);
xh(:,1) = x_Initial;
yh(:,1) = C*xh(:,1)+D*u(:,1);
e(:,1) = y(:,1)-yh(:,1);
% ▃▃ 模擬過程 ▃▃
for i = 2:Num_Sample
    if Predict_OKID == 0 % 非預測型
        xh(:,i) = G*xh(:,i-1)+H*u(:,i-1)-Lo*e(:,i-1);
    elseif Predict_OKID == 1 % 預測型
        xh(:,i) = Gd*xh(:,i-1)+Hd*u(:,i-1)+Lo*y(:,i);
    end
    yh(:,i) = C*xh(:,i)+D*u(:,i);
    e(:,i) = y(:,i)-yh(:,i);
    th(:,i) = (i-1)*Ts;
end
% ▃▃ 比對狀態維度及筆數 ▃▃
x_in = xh(:,10);
disp('   Row Col')
disp(['xh',char(9),num2str(size(xh)),char(10),...
      'ue',char(9),num2str(size(u)),char(10),...
      'yh',char(9),num2str(size(yh)),char(10),...
      'th',char(9),num2str(size(th)),char(10)])