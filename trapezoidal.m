%Trapezoidal
%Inputs

function [X Y] = trapezoidal(t,A,B,C,U)



dt=t(2)-t(1);
t0=t(1);
tf=t(end);
[N_Nodes,~]=size(A);


X(:,1)=zeros(N_Nodes,1);
for iter=1:length(t)-1
    Bstar=(1/2)*dt*(B*U(iter+1)+A*X(:,iter)+B*U(iter))+X(:,iter);
    X(:,iter+1)=(eye(N_Nodes)-dt/2.*A)\Bstar;
    Yout=C*X(:,iter);
    Y(:,iter)=Yout;
end

Y=[0 Y];
