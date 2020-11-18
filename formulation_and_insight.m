% % this script is to record the formulation in the NESOI paper
% % and illustrate the deriavation and data-shape of every variables;
% % this script could debug because of lacking key parameter for NESOI model
%% parameter
n=70;       %channels of eeg 
d=8196;     %cortex mesh for the surface of generated cortex
s=121*881;  %time course of eeg
k=20;       %the number of indepent component 
%% fomulation and every shape of variabel,note that it maybe not need in practice algorithm
%% theta1=zeros(n,s);
e1=randn(n,1);
T=randn(1,s);
theta1=kron(e1,T);
%% theta2=zeros(d,s);
e2=randn(d,1);
T2=randn(1,s);
theta2=kron(e2,T2);
%% model shape;
L=zeros(n,d);      %gain matrix
Y=zeros(n,s);      %EEG signal
% noise covariant 
C1=zeros(n,n);     
C2=zeros(d,d);
 
%% first level
Y=L*theta+theta1;              %(1)
%e1~N(0,C1);
%% second level
theta=0+theta2;                %(2)  
% e2~N(0,C2);
%% hyper para
gammma=zeros(k,1);  %hyperparameter we want to obtain to decide which fMRI's sIC contributes most.
V=cell(1,k);        %V=[V1 V2 ... Vk] ,in which Vi=zeros(d,d)
alpha=randn(1);
%% noise relationship
% C1
C1=alpha^-1*eye(n);
% C2
for i=1:k
    C2=C2+gamma(i)*V{i};
end

