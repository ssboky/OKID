%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ver. �G2013092612
% Desc.�G�Ѿ\�Ӹ�Ƨ� OKID.pdf�C
% Info.�GID_u,ID_y ����X�J Data, q �� Markov �ѼơC
%        Mult_Num = �\ = �], D_exist �� D �������L�C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G,H,C,D,Lo,Singu,Vn,Sn] = Auxi_OKID_JXL(ID_u,ID_y,q,Mult_Num,D_exist)
% Singu_Value = diag(Singu); % ��K�ƻs�ϥ�(���X�_����)
% �d�d ��X�J���׽T�w �d�d
[m,u_l] = size(ID_u); % m ����J�Ӽ�
[p,y_l] = size(ID_y); % p ����X�Ӽ�
% �d�d �����Ѽƹw�]�ȳ]�w(�q�{�]�w) �d�d
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
% �d�d �N��X�J ID_Data �H�˩�m���覡�զ� 2.7e �� �d�d
v_bar = [ID_u; ID_y]; % 2.4 ��
V_bar = [];
for i = 1:q
    V_bar = [v_bar(:,i:end-((q+1)-i)); V_bar]; % 2.7e ��
end
% �d�d D �������L�A�M�w 2.7e ���� u(q),... ���@�C���h�d �d�d
% �d�d �Ѧ� 2.7e ���P 2.34b �������t���� �d�d
y_bar = ID_y(:,q+1:end); % 2.7b ��
if (D_exist == 0)
    Y_bar = y_bar*pinv(V_bar);
    D = zeros(p,m);
    Y_bar = [D,Y_bar];
else
    V_bar = [ID_u(:,q+1:end); V_bar]; % �W�[ u(q),... ���@�C
    Y_bar = y_bar*pinv(V_bar);
    D = Y_bar(:,1:m);    
end
Y_bar(:,1:m) = []; % �D�X D ����R���A��l�� 2.10b
% �d�d �N Y_bar ���Ѭ� Y_bar_1 �P Y_bar_2 �d�d
for i = 1:q
    Y_bar_1(:,:,i) = Y_bar(:,(i-1)*(p+m)+1:(i-1)*(p+m)+m);
    Y_bar_2(:,:,i) = -Y_bar(:,(i-1)*(p+m)+(m+1):(i-1)*(p+m)+(m+p));    
end
% �d�d 2.12 �� Yk ���p�� �d�d
% �d�d ���F�K��p��A�N Yk �� k �]�� k+1�A�䤤�O Y1 �� D �d�d
% �d�d ����A�N k+1 �]�^ k (�{�� Y(:,:,1) = []) �d�d
Y(:,:,1) = D;
for k = 1:2*(Mult_Num+1)
    Left_temp = []; Right_temp = []; Sum_temp = [];
    if (k <= q) % 2.12b ��
        for j = 1:k
            Left_temp = [Left_temp Y_bar_2(:,:,j)];
            Right_temp = [Y(:,:,j); Right_temp];
        end
        Sum_temp = Left_temp*Right_temp;
        Y(:,:,k+1) = Y_bar_1(:,:,k)-Sum_temp;
    else % 2.12c ��
        for j = 1:q
            Left_temp = [Left_temp Y_bar_2(:,:,j)];
            Right_temp = [Right_temp; Y(:,:,k-j+1)];
        end
        Sum_temp = Left_temp*Right_temp;
        Y(:,:,k+1) = -Sum_temp;
    end
end
Y(:,:,1) = []; % �R�� D �Ӷ�
% �d�d 2.14 �� Yok ���p�� �d�d
Y_o(:,:,1) = Y_bar_2(:,:,1);
for k = 2:2*(Mult_Num+1)
    Left_temp = []; Right_temp = []; Sum_temp = [];
    if (k <= q) % 2.14b ��
        for j = 1:(k-1)
            Left_temp = [Left_temp Y_bar_2(:,:,j)];
            Right_temp = [Y_o(:,:,j); Right_temp];
        end
        Sum_temp = Left_temp*Right_temp;
        Y_o(:,:,k) = Y_bar_2(:,:,k)-Sum_temp;
    else % 2.14c ��
        for j = 1:q
            Left_temp = [Left_temp Y_bar_2(:,:,j)];
            Right_temp = [Right_temp; Y_o(:,:,k-j)];
        end
        Sum_temp = Left_temp*Right_temp;
        Y_o(:,:,k) = -Sum_temp;
    end
end
% �d�d 2.15 �����p�� �d�d
H_bar = [];
for i = 1:Mult_Num+1
    H_bar_temp = [];    
    for j = 1:Mult_Num+2
        temp = [Y(:,:,(i-1)+j) Y_o(:,:,(i-1)+j)];
        H_bar_temp = [H_bar_temp temp];
    end
    H_bar = [H_bar; H_bar_temp];
end
H_bar_0 = H_bar(:,1:(p+m)*(Mult_Num+1)); % �� k = 1 �ɡAH_bar ���d��
H_bar_1 = H_bar(:,(p+m)+1:end); % �� k = 2 �ɡAH_bar ���d��
% �d�d H_bar_0 ���_���Ȥ��� �d�d
[V,Singu,S] = svd(H_bar_0);
n_min = q*p;
Vn = V(:,1:n_min);
Sn = S(:,1:n_min);
Singun = [Singu(1:n_min,1:n_min)];
Singu_Value = diag(Singu);
% �d�d 2.17a �� �d�d
G = (Singun^-0.5)*Vn'*H_bar_1*Sn*(Singun^-0.5); 
% �d�d 2.17b �� �d�d
InputMatrix_temp = (Singun^0.5)*Sn';
H = InputMatrix_temp(:,1:m); 
Lo = InputMatrix_temp(:,(m+1):(m+p));
% �d�d 2.17c �� �d�d
OutputMatrix_temp = Vn*(Singun^0.5);
C = OutputMatrix_temp(1:p,:);