function [xout sysout q Matrices] = KrylovMOR(K, M, C, Forces, ObservationVector, Time, U, X0, Order )


%References: 

%Krylov subspace techniques for reduced-order modeling of large-scale
%dynamical systems
%Zhaojun Bai
%Works for one observation only
%see also: doi:10.1088/1757-899X/10/1/012118


[n_dof,~]=size(K);

q=zeros(n_dof,Order);
p=zeros(n_dof,Order);

r0=Forces\K;
r0=r0/norm(r0);
q(:,1)=r0;
p(:,1)=zeros(n_dof,1);
for j=1:Order
    r=K\C*q(:,j)+K\M*p(:,j);
    s=q(:,j);
    for i=1:j
        t(i*j)=q(:,i)'*r;
        r=r-q(:,j)*t(i*j);
        s=s-p(:,j)*t(i*j);
    end
    t(j+j)=norm(r,2);
    if t(j+j)==0
        break
    else
        q(:,j+1)=r/t(j+j);
        p(:,j+1)=s/t(j+j);
    end
end

Mr=q'*M*q;
Cr=q'*C*q;
Kr=q'*K*q;
Fr=q'*Forces;


Xd0=0;
Xdd0=0;

xr=q'*X0;
xdr=q'*Xd0;
xddr=q'*Xdd0;


Xr=NewmarkIntegrate(Kr,Mr,Cr,Fr,Time,U,xr,xdr,xddr,false);

c=q'*ObservationVector;
xout=c'*Xr;


A=[zeros(size(Kr)),eye(size(Kr));...
    -Mr\Kr, -Mr\Cr];

B=[zeros(Order+1,1);Mr\-Fr];

C=[zeros(size(c)); c]';

sysout=ss(A,B,C,0);

Matrices=struct;
Matrices.Mr=Mr;
Matrices.Cr=Cr;
Matrices.Kr=Kr;
Matrices.Fr=Fr;
