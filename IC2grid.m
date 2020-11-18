clear all; clc;
%% 一些说明
% 顶点：cortex mesh ，指的是大脑皮质的表面顶点（surface mesh）
% 体素：voxel ，指的是fMRI信号每个最小的单元
% 脑体积坐标：volume coordinate，指的是fMRI按照一定的标号顺序来表征体素在空间的位置
% sIC：空间ICA分解。得到的是独立成分下的体素激活信息。同时做过稀疏化处理，因此有
% brain_ind：指的是sIC的列向量中每个序号对应脑体积坐标的哪个序号。
% Vox2mesh:体素与皮质网格的传递矩阵
% Q:高斯模糊后的邻接矩阵
%% green matrix and adjancent matrix
load(['./data/','source_inf.mat']);
Green=cort.VertConn;
d=size(Green,1);
Q=zeros(d,d);
rho=0.6;
for i=1:8
    Q=Q+(rho/factorial(i)*(Green^i));
end
save Qmatrix.mat Q
%% co-register fMRI and EEG
load(['./data/','fmri-ica.mat']);
% ica=compSet;
% clear compSet
sIC=ic;
load(['./data/','mask_ind']);
brain_ind=mask_ind;                       %sIC对应的体素索引
clear mask_ind
%[X Y Z]=ind2sub([79,95,79],brain_id);   %得到fMRI的volume坐标
%% 数据形状信息
k=size(sIC,1);
d=size(Green,1);    
%% z-scores 
u=mean(sIC,2);              %均值
sigma=sqrt(var(sIC,0,2));   %等价于std。方差
sIC_z=(sIC-u)./sigma;       %Z-transform
zscore=abs(sIC_z);          %负值也考虑进去？
clear u sigma sIC_z
%% 确定sIC下满足条件的体素位置索引
sIC_ind=cell(k,1);       %确定sIC下满足条件的体素位置索引
for i=1:k
    sIC_ind{i}=find(zscore(i,:)>3);   %index用于存储每个独立成分下z分数大于3的voxel值
end
%co-register voxel and get U matrix
vox2mesh=cort.tess2mri_interp;      %体素与皮质网格的传递矩阵
[Y mesh_ind]=max(vox2mesh,[],1);        %每个顶点对应的最大可能性的体素位置
clear Y
%mesh_ind 表示的是每个皮质顶点对应的最有可能的体素值
%% 将满足条件的sIC的位置与体素对准
sIC2vox=cell(k,1);                       
for i=1:k;
    sIC2vox{i}=brain_ind(sIC_ind{i},1);          %找出满足zscore>3的sIC所对应的体素位置
end
%% 将sIC的体素与皮质顶点建立关系： mesh_ind <--->sIC2vol
% plan 1 :用哈希表
%通过哈希表构建索引与数值的关系：
% volume_ht=java.util.Hashtable;
% for i=1:d
%     volume_ht.put(mesh_ind(1,i),i);
% end
% 初始化索引矩阵
non0index=cell(k,1);
% mesh_list=cell(k,1);
ICAgrid=zeros(d,k);
% 计算ICA格电矩阵
for j=1:k
% m=size(sIC2vox{j},1);
% mesh_list{j}=zeros(m,1);
% for i=1:m
%     %% 判断为空，输入0值。
%     if isempty(volume_ht.get(sIC2vox{j}(i)))
%         mesh_list{j}(i)=0;
%     else
%         mesh_list{j}(i)=volume_ht.get(sIC2vox{j}(i));
%     end
%     %mesh list指的是每个sIC对应的顶点坐标。
%     
%     %%  再从中选出不为0的坐标，意味着这些sIC的体素点是能够配准上顶点的
%     non0index{j}=mesh_list{j}(find(mesh_list{j}~=0),1);
% end
   %plan 2 直接用筛选函数一句话搞定
    [non0index{j} ia ib]=intersect(mesh_ind',sIC2vox{j});
    ICAgrid(ia,j)=1;   %按照索引将激活的顶点赋值为1
end
sparse(ICAgrid);
ic_num=sum(ans,1)
save ICA_grid ICAgrid
clear j m i
