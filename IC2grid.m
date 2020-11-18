clear all; clc;
%% һЩ˵��
% ���㣺cortex mesh ��ָ���Ǵ���Ƥ�ʵı��涥�㣨surface mesh��
% ���أ�voxel ��ָ����fMRI�ź�ÿ����С�ĵ�Ԫ
% ��������꣺volume coordinate��ָ����fMRI����һ���ı��˳�������������ڿռ��λ��
% sIC���ռ�ICA�ֽ⡣�õ����Ƕ����ɷ��µ����ؼ�����Ϣ��ͬʱ����ϡ�軯���������
% brain_ind��ָ����sIC����������ÿ����Ŷ�Ӧ�����������ĸ���š�
% Vox2mesh:������Ƥ������Ĵ��ݾ���
% Q:��˹ģ������ڽӾ���
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
brain_ind=mask_ind;                       %sIC��Ӧ����������
clear mask_ind
%[X Y Z]=ind2sub([79,95,79],brain_id);   %�õ�fMRI��volume����
%% ������״��Ϣ
k=size(sIC,1);
d=size(Green,1);    
%% z-scores 
u=mean(sIC,2);              %��ֵ
sigma=sqrt(var(sIC,0,2));   %�ȼ���std������
sIC_z=(sIC-u)./sigma;       %Z-transform
zscore=abs(sIC_z);          %��ֵҲ���ǽ�ȥ��
clear u sigma sIC_z
%% ȷ��sIC����������������λ������
sIC_ind=cell(k,1);       %ȷ��sIC����������������λ������
for i=1:k
    sIC_ind{i}=find(zscore(i,:)>3);   %index���ڴ洢ÿ�������ɷ���z��������3��voxelֵ
end
%co-register voxel and get U matrix
vox2mesh=cort.tess2mri_interp;      %������Ƥ������Ĵ��ݾ���
[Y mesh_ind]=max(vox2mesh,[],1);        %ÿ�������Ӧ���������Ե�����λ��
clear Y
%mesh_ind ��ʾ����ÿ��Ƥ�ʶ����Ӧ�����п��ܵ�����ֵ
%% ������������sIC��λ�������ض�׼
sIC2vox=cell(k,1);                       
for i=1:k;
    sIC2vox{i}=brain_ind(sIC_ind{i},1);          %�ҳ�����zscore>3��sIC����Ӧ������λ��
end
%% ��sIC��������Ƥ�ʶ��㽨����ϵ�� mesh_ind <--->sIC2vol
% plan 1 :�ù�ϣ��
%ͨ����ϣ������������ֵ�Ĺ�ϵ��
% volume_ht=java.util.Hashtable;
% for i=1:d
%     volume_ht.put(mesh_ind(1,i),i);
% end
% ��ʼ����������
non0index=cell(k,1);
% mesh_list=cell(k,1);
ICAgrid=zeros(d,k);
% ����ICA������
for j=1:k
% m=size(sIC2vox{j},1);
% mesh_list{j}=zeros(m,1);
% for i=1:m
%     %% �ж�Ϊ�գ�����0ֵ��
%     if isempty(volume_ht.get(sIC2vox{j}(i)))
%         mesh_list{j}(i)=0;
%     else
%         mesh_list{j}(i)=volume_ht.get(sIC2vox{j}(i));
%     end
%     %mesh listָ����ÿ��sIC��Ӧ�Ķ������ꡣ
%     
%     %%  �ٴ���ѡ����Ϊ0�����꣬��ζ����ЩsIC�����ص����ܹ���׼�϶����
%     non0index{j}=mesh_list{j}(find(mesh_list{j}~=0),1);
% end
   %plan 2 ֱ����ɸѡ����һ�仰�㶨
    [non0index{j} ia ib]=intersect(mesh_ind',sIC2vox{j});
    ICAgrid(ia,j)=1;   %��������������Ķ��㸳ֵΪ1
end
sparse(ICAgrid);
ic_num=sum(ans,1)
save ICA_grid ICAgrid
clear j m i
