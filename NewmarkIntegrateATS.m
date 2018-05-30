function X = NewmarkIntegrate(K,M,C,F,t,U,X0,Xd0,Xdd0,animate,animatestruct)


dt=t(2)-t(1);
[n_dof,~]=size(K);

del=0.5;
a=0.25*(0.5+del)^2;
a0=1/(a*dt^2);
a1=del/(a*dt);
a2=1/(a*dt);
a3=1/(2*a)-1;
a4=del/a-1;
a5=(dt/2)*(del/a-2);
a6=dt*(1-del);
a7=del*dt;

if nargin>10
    Nodes=animatestruct.nodes;
    Elements=animatestruct.elements;
    DOF=animatestruct.dof;
end


Kef=K+a0*M+a1*C;
%[L,U,P]=lu(Kef,'vector');

X=X0;
Xd=Xd0;
Xdd=Xdd0;
Reff=zeros(n_dof,1);
framecounter=1;
while thistime<t(end);
    Reff(:,i+1)=U(i+1)*F+M*(a0*X(:,i)+a2*Xd(:,i)+a3*Xdd(:,i))+C*(a1*X(:,i)+a4*Xd(:,i)+a5*Xdd(:,i));
    Xnew=Kef\Reff(:,i+1);
    X(:,i+1)=Xnew;
    Xdd(:,i+1)=a0*(X(:,i+1)-X(:,i))-a2*Xd(:,i)-a3*Xdd(:,i);
    Xd(:,i+1)=Xd(:,i)+a6*Xdd(:,i)+a7*Xdd(:,i+1);
    
    if animate==true
        if framecounter==50
            framecounter=1;
        [thisframeNodes,thisframeCS_Origin]=newNodes(Nodes,Xnew,Elements,DOF);        
        PlotInitial(thisframeNodes,Elements,thisframeCS_Origin,false);
        titlestr=sprintf('Transient Response t = : %8.3f', t(i));

        title(titlestr);
        
        drawnow;
        end
    end
   framecounter=framecounter+1 ;
end



