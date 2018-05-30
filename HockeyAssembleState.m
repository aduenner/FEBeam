function [A,B,C] = HockeyAssembleState(K,M,C,Forces)

n_dof=length(Forces);
%


A=[zeros(size(K)),eye(size(K));...
    -M\K, -M\C];

B=[zeros(n_dof,1);M\-Forces];

C=[zeros(n_dof,n_dof),diag(ones(1,n_dof))];




