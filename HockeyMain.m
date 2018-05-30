%Generate Points
clc
clear all
close all
set(0, 'DefaultFigureRenderer', 'OpenGL');
set(0,'DefaultFigureWindowStyle','docked')

%%%%%%%%%%%%%%%%%%%%%%%%  USER OPTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Geometry
    wristStick = false;
%Input
      inputType =3; %1--Gaussian  2--Step  3--Impulse
 initialPreload = false;
%Damping Params
           freq1 = 3;
           freq2 = 25;
        damping1 = 0.05;  
        damping2 = 0.1;
 
%Nodes
      ShaftNodes = 20;
      BladeNodes = 20;        
%Force (multiplied by the gaussian distribution of U later on)
       forcenode = ShaftNodes+BladeNodes;
    forcenodedir = 3;
   forceConstant = 7;
 
%Solve Options
updateStiffness = false;
  explicitsolve = false;
      Integrate = true;

%Animate Options
    liveanimate = true;
    animatepost = false;

%Frequency Domain Display Options
       showBode = false;
      showEigen = true;
EigenModeNumber = 5;
    staticEigen = true;

%Model Reduction
 CondensedModel = false;
          Order = 2;
    reduceModel = true;


%Time
              t0 = 0;
              dt = 1e-4;
              tf = 1;
               t = t0:dt:tf;



%%%%%%%%%%%%%%%%%%%%%%%%%%%  END OF OPTIONS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get Initial Nodes, Elements, Element Properties given n1 nodes in region
%of shaft and n2 nodes in region of blade

[alpha, beta]=rayleighCoefficients(freq1,damping1,freq2,damping2);

if wristStick==true
    [Nodes,Nodal_DOFS,Elements,El_Properties]=WristHockeyInitialNodes(ShaftNodes,BladeNodes);
    forceConstant=forceConstant*2;
else
    [Nodes,Nodal_DOFS,Elements,El_Properties]=HockeyInitialNodes(ShaftNodes,BladeNodes);
end

%Get Initial Orientation
Element_GCS=getElementOrientation(Nodes,Elements);

%Initial K,M,C global matrices
[DOF_Legend,Kg,Mg,Cg]=AssembleMatrices(Nodal_DOFS,Elements,El_Properties,Element_GCS,alpha,beta);
[n_dof,~]=size(Kg);
Forces=zeros(n_dof,1);
ForceIndexSearch=[forcenode forcenodedir];
[~,ForceNodeIndex]=ismember(ForceIndexSearch,DOF_Legend(:,2:3),'rows');

Forces(ForceNodeIndex)=forceConstant;

%State Space represntation
[A,B,C] = HockeyAssembleState(Kg,Mg,Cg,Forces);

switch inputType
    case 1 %gaussian
        U=zeros(size(t));
        U(1)=0;
        tistart=5;
        tiend=15;
        impacttime=t(tistart:tiend);
        impactcenter=t(tiend-tistart+1);
        Uimpact=normpdf(impacttime,impactcenter,(t(tiend)-t(tistart))/4);
        Uimpact=Uimpact/max(Uimpact)*80;
        U(2:length(impacttime)+1)=Uimpact;
        U=U';
    
    case 2 %Step Input
        U=ones(size(t));
        U(1)=0;
    
    case 3 %Impulse
        U=zeros(size(t));
        U(2)=1000;
        
    otherwise
        U=zeros(size(t));
end

%Initialize Output

if initialPreload==true
    tic
    X0=Kg\Forces;
    timesimstatic=toc
    [Nodes,Element_GCS]=newNodes(Nodes,X0,Elements,DOF_Legend);
else
    X0=zeros(n_dof,1);
end
Xd0=X0;
Xdd0=X0;

if reduceModel==true
    ObservationVector=zeros(size(Forces));
    ObservationVector(ForceNodeIndex)=1;
        tic
    [xout redsysout q reducedMatrices] = KrylovMOR(Kg, Mg, Cg, Forces, ObservationVector, t, U, X0, Order );
        timesimreducedsys=toc
    if CondensedModel==true

    for i=1:length(qsum)
        killnode=DOF_Legend(i,2);
        killdofindex=DOF_Legend(i,3)+1;
         if abs(qsum(i))<1e-6
            Nodal_DOFS(killnode,killdofindex)=0;
        else
        end
    end
    
    [DOF_Legend,Kg,Mg,Cg]=AssembleMatrices(Nodal_DOFS,Elements,El_Properties,Element_GCS,alpha,beta);
    [n_dof,~]=size(Kg);
    Forces=zeros(n_dof,1);
    [~,ForceNodeIndex]=ismember(ForceIndexSearch,DOF_Legend(:,2:3),'rows');
    Forces(ForceNodeIndex)=forceConstant;
    
        if initialPreload==true
            X0=Kg\Forces;
            [Nodes,Element_GCS]=newNodes(Nodes,X0,Elements,DOF_Legend);
        else
            X0=zeros(n_dof,1);
        end
        Xd0=zeros(size(X0));
        Xdd0=Xd0;
    end
end

if explicitsolve==true
    A0=(1/dt^2*Mg+1/dt*Cg);
    A1=(2/dt^2*Mg-Kg);
    A2=-1/dt^2*Mg+1/(2*dt)*Cg;
    X=X0;
    Xd=X0;
    Xdd=X0;
    X(:,2)=X0-dt*Xd(:,1)+dt^2/2*Xdd(:,1);
    for i=2:length(t)-1
        Bstar=A1*X(:,i)+A2*X(:,i-1)+U(i).*Forces;
        Xnew=A0\Bstar;
        X(:,i+1)=Xnew;
        
        if updateStiffness==true
            [postNodes,postElement_CS]=newNodes(Nodes,X0,Elements,DOF_Legend);
            [~,Kg,Mg,Cg]=AssembleMatrices(Nodal_DOFS,Elements,El_Properties,postElement_CS,alpha,beta);
            A0=(1/dt^2*Mg+1/dt*Cg);
            A1=(2/dt^2*Mg-Kg);
            A2=-1/dt^2*Mg+1/(2*dt)*Cg;
        else
        end
        
    end
    
elseif Integrate==true
    animatestruct=struct;
    animatestruct.nodes=Nodes;
    animatestruct.elements=Elements;
    animatestruct.dof=DOF_Legend;
    if liveanimate==true
        
        X = NewmarkIntegrate(Kg,Mg,Cg,Forces,t,U,X0,Xd0,Xdd0,true,animatestruct);
        pause
    else
        tic
        X = NewmarkIntegrate(Kg,Mg,Cg,Forces,t,U,X0,Xd0,Xdd0,false,animatestruct);
        timesimfullmodel=toc
    end
end



if animatepost==true
    figure()
    videoout=VideoWriter('HockeyStickGaussianOnlyHighDamping');
    videoout.FrameRate=30;
    open(videoout);
    for i=1:length(X)/50
        j=i*50;
        nodedisp=X(:,j);
        [thisframeNodes,thisframeCS_Origin]=newNodes(Nodes,nodedisp,Elements,DOF_Legend);
        PlotInitial(thisframeNodes,Elements,thisframeCS_Origin,false);
        drawnow;
        mframe=getframe;
        writeVideo(videoout,mframe);
        
    end
    close(videoout);
else
end

if showEigen==true
    figure()
    [V,D]=eig(A);
    for i=1:5
        EigenModeNumber=2*(i-1)+1;
        [EV,EVindex]=HockeyEigen(A,Nodes,Elements,DOF_Legend,EigenModeNumber);
        EV=EV./(2*pi);
        pause
    end
end

if staticEigen==true
    figure()
    for i=1:4
        subplot(2,2,i)
        Vstatic=V(end/2+1:end,EVindex(2*i-1));
        Vstatic=real(Vstatic)/max(abs(real(Vstatic)));
        Vstatic=Vstatic/max(abs(Vstatic(3:3:end)))*0.2;
        [n1,f1]=newNodes(Nodes,Vstatic,Elements,DOF_Legend);
        PlotInitial(n1,Elements,f1,false,[-0.3 1 -2 0.1 -0.3 0.3]);
        hold on
        PlotInitial(Nodes,Elements,Element_GCS,false,[-0.3 1 -2 0.1 -0.3 0.3]);
        title(['Eigenmode ',num2str(i),' \lambda=',num2str(EV(2*i-1)),'hz']);
        pause
    end
end


if Integrate==true || explicitsolve==true
    yend=X(end-3,:);
    %plot(t,yend);
else
end

if reduceModel==true
    figure()
    Cbode=zeros(1,2*n_dof);
    Cbode(end-3)=1;
    sysout=ss(A,B,Cbode,0);
    bodefreqs={0.1,1000};
    bode(sysout,bodefreqs)
    hold on
    title('Blade Frequency Response');
    h = gcr;
    setoptions(h,'FreqUnits','Hz')
    setoptions(h,'PhaseMatching','on')
    set(gcf,'color','w')

        bode(redsysout,bodefreqs)
        legend('Reduced Model','Original')
    pause
end

if Integrate==true && reduceModel==true
    figure()
    plot(t,yend)
    hold on
    plot(t,xout)
    legend('Reduced Model','Original')
    title('Blade Z Displacement vs Time');
    xlabel('time')
    ylabel('Blade Z Displacement (m)');
    set(gcf,'color','w');
    rmsd=sqrt(sum((xout-yend).^2)/length(yend))
    pause
end



