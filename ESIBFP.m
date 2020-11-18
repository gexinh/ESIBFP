%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&&&&&&
%%%%%%%%%%%%%%%   the algorithm   &&&&&&&&&&&&&&&&&&&&&&
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&&&&&&&&
%% load data
load('./data/eeg-ica.mat');
Y=topography;
load('./data/gainm.mat');
L=gainmatrix;
load('ICA_grid');
U=ICAgrid;
load('Qmatrix.mat')
Q=Q;                 %just in order to illustrate where the Q matrix comes from

%% select which plan to caculate
 p=3;
 load(['covmatrix',num2str(p)]);

% load('./data/eeg_select')
% Y=Y(:,eegselect);
%% fetch the data size
% n:channels
% s;the number of tIC
% d;cortex 
% k:the number of sIC
[n s]=size(Y);
d=size(L,2);
k=size(U,2);
%% initial the parameters
% gamma=zeros(k,1);
% alpha=rand(1);
% C1=eyes(n); 
% C2=zeros(d,d);             
%% V matrix=sum(U*q*q')
%%%%%%%%%%%%%%%%%%%%%%%%%%% plan 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%
% V1=eye(n);
% V={};
% for j=1:k
%     index=find(U(:,j)==1);
% %     Q()
%     % plan 1
%     q=sum(Q(:,index),2)./sqrt(length(index));
%     V{j}=diag(q);
%     %�����������ͳ��ÿ�����Լ����������ƽ������ֵ��   
% end
% for i = 1:length(V)      %lenth=sIC=20
%     LVL{i} = L*V{i}*L';  %L*V*L'
% end
% Q  = [V1 LVL];
% save covmatrix1.mat Q
% load covmatrix1

%%%%%%%%%%%%%%%%%%%%%%%%% plan 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V1=eye(n);
% V={};
% for j=1:k
%     index=find(U(:,j)==1);
% %     % plan 2
%      q=Q(:,index)*Q(:,index)';
%      V{j}=q./length(index);
%     
% end
% 
% % ��������ֵ
% for i = 1:length(V)      %lenth=sIC=20
%     LVL{i} = L*V{i}*L';  %L*V*L'
% end
% Q  = [V1 LVL];
% save covmatrix2.mat Q
% 
% load covmatrix2

%%%%%%%%%%%%%%%%%%%%%% plan 3 for cherish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     V1=eye(n);
%     V={};
%     LVL = cell(0);
%     sum_C2d = 0;
%     for i = 1:k
%         % GM ���ڽӾ���
%         % fmriic����vertice������sICA�ɷ�
%         index=find(U(:,i)==1);
%         Q_temp = Q(:,index);
%         Q_temp = max(Q_temp');  %�ҵ��óɷ�������һ������
%         Q_temp = Q_temp';
%         diag_C2 = Q_temp / sum(Q_temp);
%        
%         actN = d - sum(diag_C2 == 0) ;  %�ҳ��������ܸ���
%         fprintf(' %d    :  %d  \n', i, actN);
%       %  sum_C2d = sum_C2d + diag_C2;        %����Ϊ����MSP������ͳ��
%         V{i} = diag(diag_C2); 
%         LVL{i} = sparse(L * V{i} * L');       
%         
%     end
% Q  = [V1 LVL];
% save covmatrix3.mat Q
% 
% load covmatrix3

% % ���ƻ�������

%% set model to train
hp_matrix=zeros(s,k+1);
for i=1:s
    %% process Y signal
    Y_mean=mean(Y,1);
    Y_ct=Y(:,i)-Y_mean(:,i);
    YY=Y_ct(:,1)*Y_ct(:,1)';
    %% process covariance cell
 
    [Cy,h,ph,F]=spm_reml_sc(YY,[],Q);  %�źŷ���������
    hp_matrix(i,:)=h;  %�Զ�ת��
%CyΪgamma*V����ΪЭ�������
%hΪgamma����
%phΪgamma�ľ��Ⱦ���
%FΪ������
%Q����ƺõľ���
%Y��ԭʼ�ź�
end
%% get trained parameters and plot
% C2=zeros(m);
% for i = 1:length(Qf)
%     C2 = C2+h(length(Qt)+i)*Qf{i};
% end
% 
% C1=zeros(n);
% for i = 1:length(Qt)
%     C1 = C1+h(i)*Qt{i};
% end
% 
% % w=inv(G'*pinv(C1)*G+pinv(C2))*G'*pinv(C1)*Y;
% w=C2*G'*inv(G*C2*G'+C1)*Y;
% h=h(length(Qt)+1:end);

hp_matrix = hp_matrix ./ max(hp_matrix,[],2);
h=hp_matrix(:,2:k + 1);
h(:,[3 2])=h(:,[2 3]);
save hmatrix h;
plotConfMat(h);
% surf(hp_matrix(:,2:k + 1));                          %չʾЧ��