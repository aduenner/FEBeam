function X = NewmarkIntegrate(Matrices,dt,q)

exitnow=false;
U=0
function crazy_gui
f = figure('KeyPressFcn', @(handle,event) disp(event.Key), 'Position',[500,500,500,500]);
hButton = uicontrol('Style','pushbutton','String','Crazy Button',...
        'Position',[200,200,100,100],...
        'Callback',@button_Callback);
     h = animatedline
function button_Callback(handle,event)      
  U=U+1
end
end



M=Matrices.Mr;
C=Matrices.Cr;
K=Matrices.Kr;
F=Matrices.Fr;


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


Kef=K+a0*M+a1*C;
%[L,U,P]=lu(Kef,'vector');
X0=zeros(n_dof,1);

X=X0;
Xd=X0;
Xdd=X0;
Reff=zeros(n_dof,1);

          
i=1;
 % detect mouse click
while exitnow==false
    
    Reff(:,i+1)=U*F+M*(a0*X(:,i)+a2*Xd(:,i)+a3*Xdd(:,i))+C*(a1*X(:,i)+a4*Xd(:,i)+a5*Xdd(:,i));
    Xnew=Kef\Reff(:,i+1);
    X(:,i+1)=Xnew;
    Xdd(:,i+1)=a0*(X(:,i+1)-X(:,i))-a2*Xd(:,i)-a3*Xdd(:,i);
    Xd(:,i+1)=Xd(:,i)+a6*Xdd(:,i)+a7*Xdd(:,i+1);
    

     addpoints(h,i*dt,F'*Xnew);
    drawnow
           i=i+1
    end

end



