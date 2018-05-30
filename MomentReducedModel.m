function [A_,B_,C_,Vq] = MomentReducedModel(A,B,C,order)

tbr=true;

Vq(:,1)=A\B;
    Vqnorms(1)=norm(Vq(:,1),2);
    Vq(:,1)=Vq(:,1)/Vqnorms(1);
    for qiter=2:order
        Vq(:,qiter)=A\Vq(:,qiter-1);
        for ortholoop=1:qiter-1
            beta = Vq(:,qiter)'*Vq(:,qiter-ortholoop);
            Vq(:,qiter)=Vq(:,qiter)-beta*Vq(:,qiter-ortholoop);
        end
        Vqnorms(qiter)=norm(Vq(:,qiter),2);
        Vq(:,qiter)=Vq(:,qiter)/Vqnorms(qiter);

    end
        
        sysout=ss(A,B,C,0);
       sysreduced=balancmr(sysout,order);
    if tbr==true
    A_= sysreduced.A;
    B_=sysreduced.B;
    C_=sysreduced.C;
    else
    A_= Vq'*A*Vq;
    B_=Vq'*B;
    C_=Vq'*C';
    C_=C';
    end
%     
