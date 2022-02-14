%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ver. ：2013092612
% Desc.：參閱該資料夾 OKID.pdf。
% Info.：ID_u,ID_y 為輸出入 Data, q 為 Markov 參數。
%        Mult_Num = α = β, D_exist 為 D 項的有無。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G,H,C,D,Lo,Singu,Vn,Sn] = Auxi_OKID_JXL(ID_u,ID_y,q,Mult_Num,D_exist)
% Singu_Value = diag(Singu); % 方便複製使用(取出奇異值)
% ▃▃ 輸出入維度確定 ▃▃
[m,u_l] = size(ID_u); % m 為輸入個數
[p,y_l] = size(ID_y); % p 為輸出個數
% ▃▃ 相關參數預設值設定(默認設定) ▃▃
Check_q = exist('q');
if (Check_q == 0 | q == 0)
    q = 1;
end
% -----
Check_Mult_Num = exist('Mult_Num');
if (Check_Mult_Num == 0 | Mult_Num == 0)
    Mult_Num = 2;
end
% -----
Check_D_exist = exist('D_exist');
if (Check_D_exist == 0 | D_exist == 0)
    D_exist = 0;
end
% ▃▃ 將輸出入 ID_Data 以倒放置的方式組成 2.7e 式 ▃▃
v_bar = [ID_u; ID_y]; % 2.4 式
V_bar = [];
for i = 1:q
    V_bar = [v_bar(:,i:end-((q+1)-i)); V_bar]; % 2.7e 式
end
% ▃▃ D 項的有無，決定 2.7e 式的 u(q),... 那一列的去留 ▃▃
% ▃▃ 參考 2.7e 式與 2.34b 式間的差異性 ▃▃
y_bar = ID_y(:,q+1:end); % 2.7b 式
if (D_exist == 0)
    Y_bar = y_bar*pinv(V_bar);
    D = zeros(p,m);
    Y_bar = [D,Y_bar];
else
    V_bar = [ID_u(:,q+1:end); V_bar]; % 增加 u(q),... 那一列
    Y_bar = y_bar*pinv(V_bar);
    D = Y_bar(:,1:m);    
end
Y_bar(:,1:m) = []; % 求出 D 項後刪除，其餘為 2.10b
% ▃▃ 將 Y_bar 分解為 Y_bar_1 與 Y_bar_2 ▃▃
for i = 1:q
    Y_bar_1(:,:,i) = Y_bar(:,(i-1)*(p+m)+1:(i-1)*(p+m)+m);
    Y_bar_2(:,:,i) = -Y_bar(:,(i-1)*(p+m)+(m+1):(i-1)*(p+m)+(m+p));    
end
% ▃▃ 2.12 式 Yk 的計算 ▃▃
% ▃▃ 為了便於計算，將 Yk 的 k 設為 k+1，其中令 Y1 為 D ▃▃
% ▃▃ 後續再將 k+1 設回 k (程式 Y(:,:,1) = []) ▃▃
Y(:,:,1) = D;
for k = 1:2*(Mult_Num+1)
    Left_temp = []; Right_temp = []; Sum_temp = [];
    if (k <= q) % 2.12b 式
        for j = 1:k
            Left_temp = [Left_temp Y_bar_2(:,:,j)];
            Right_temp = [Y(:,:,j); Right_temp];
        end
        Sum_temp = Left_temp*Right_temp;
        Y(:,:,k+1) = Y_bar_1(:,:,k)-Sum_temp;
    else % 2.12c 式
        for j = 1:q
            Left_temp = [Left_temp Y_bar_2(:,:,j)];
            Right_temp = [Right_temp; Y(:,:,k-j+1)];
        end
        Sum_temp = Left_temp*Right_temp;
        Y(:,:,k+1) = -Sum_temp;
    end
end
Y(:,:,1) = []; % 刪除 D 該項
% ▃▃ 2.14 式 Yok 的計算 ▃▃
Y_o(:,:,1) = Y_bar_2(:,:,1);
for k = 2:2*(Mult_Num+1)
    Left_temp = []; Right_temp = []; Sum_temp = [];
    if (k <= q) % 2.14b 式
        for j = 1:(k-1)
            Left_temp = [Left_temp Y_bar_2(:,:,j)];
            Right_temp = [Y_o(:,:,j); Right_temp];
        end
        Sum_temp = Left_temp*Right_temp;
        Y_o(:,:,k) = Y_bar_2(:,:,k)-Sum_temp;
    else % 2.14c 式
        for j = 1:q
            Left_temp = [Left_temp Y_bar_2(:,:,j)];
            Right_temp = [Right_temp; Y_o(:,:,k-j)];
        end
        Sum_temp = Left_temp*Right_temp;
        Y_o(:,:,k) = -Sum_temp;
    end
end
% ▃▃ 2.15 式的計算 ▃▃
H_bar = [];
for i = 1:Mult_Num+1
    H_bar_temp = [];    
    for j = 1:Mult_Num+2
        temp = [Y(:,:,(i-1)+j) Y_o(:,:,(i-1)+j)];
        H_bar_temp = [H_bar_temp temp];
    end
    H_bar = [H_bar; H_bar_temp];
end
H_bar_0 = H_bar(:,1:(p+m)*(Mult_Num+1)); % 當 k = 1 時，H_bar 的範圍
H_bar_1 = H_bar(:,(p+m)+1:end); % 當 k = 2 時，H_bar 的範圍
% ▃▃ H_bar_0 的奇異值分解 ▃▃
[V,Singu,S] = svd(H_bar_0);
n_min = q*p;
Vn = V(:,1:n_min);
Sn = S(:,1:n_min);
Singun = [Singu(1:n_min,1:n_min)];
Singu_Value = diag(Singu);
% ▃▃ 2.17a 式 ▃▃
G = (Singun^-0.5)*Vn'*H_bar_1*Sn*(Singun^-0.5); 
% ▃▃ 2.17b 式 ▃▃
InputMatrix_temp = (Singun^0.5)*Sn';
H = InputMatrix_temp(:,1:m); 
Lo = InputMatrix_temp(:,(m+1):(m+p));
% ▃▃ 2.17c 式 ▃▃
OutputMatrix_temp = Vn*(Singun^0.5);
C = OutputMatrix_temp(1:p,:);