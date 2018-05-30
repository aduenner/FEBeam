%Generate Points
clc
clear all
close all

%Damping Params
freq1=3;
freq2=25;
damping1=0.05;
damping2=0.05;
[alpha, beta]=rayleighCoefficients(freq1,damping1,freq2,damping2);
damping=.100;

%Nodes
ShaftNodes=40;
BladeNodes=40;

%Force (multiplied by the gaussian distribution of U later on)
forceConstant=100; 

%Get Initial Nodes, Elements, Element Properties given n1 nodes in region
%of shaft and n2 nodes in region of blade
[Nodes,Nodal_DOFS,Elements,El_Properties]=HockeyInitialNodes(ShaftNodes,BladeNodes);

[Ndof,~]=size(Nodes);
%Get Initial Orientation
Element_GCS=getElementOrientation(Nodes,Elements);

InitialNodes=Nodes;
Initial_GCS=Element_GCS;

%Initial K,M,C global matrices
[DOF_Legend,Kg,Mg,Cg]=AssembleMatrices(Nodal_DOFS,Elements,El_Properties,Element_GCS,alpha,beta);
[n_dof,~]=size(Kg);

%Initial displacement
%x0=Kg\Forces;
        %[thisframeNodes,thisframeCS_Origin]=newNodes(Nodes,x0,Elements,DOF_Legend);
        %PlotInitial(thisframeNodes,Elements,thisframeCS_Origin,false);

%State Space represntation
%[A,B,C] = HockeyAssembleState(Kg,Mg,Cg,Forces);
%Cbode=zeros(1,2*n_dof);
%Cbode(end-3)=1;
%sysout=ss(A,B,Cbode,0);
%[V,D]=eig(A);

%[Arows,Acols]=size(A);
%Time
t0=0;
dt=1e-5;
tf=0.5;
t=t0:dt:tf;

%Input

U=zeros(size(t));
U(1)=0;
impacttime=1e-6:1e-6:100e-6;
impactcenter=impacttime(floor(end/2));
Uimpact=normpdf(impacttime,impactcenter,5e-6)./8e4*80;
U(2:length(impacttime)+1)=Uimpact;
U=U';

%Initialize Output




%X0=Kg\Forces;
X0=gpuArray(zeros(n_dof,1));
Xd0=X0;
Xdd0=X0;

%[postNodes,postElement_CS]=newNodes(Nodes,X0,Elements,DOF_Legend);
%[~,Kg,Mg,Cg]=AssembleMatrices(Nodal_DOFS,Elements,El_Properties,postElement_CS,alpha,beta);


%Definte Forces

Forces=zeros(n_dof,1);
Forces(end-3)=forceConstant;

Forces=gpuArray(Forces);
X=X0;
Xd=X0;
Xdd=X0;
Xi=X0;
Xi1=X0;

A0=gpuArray((1/dt^2*Mg+1/dt*Cg));
A1=gpuArray((2/dt^2*Mg-Kg));
A2=gpuArray(-1/dt^2*Mg+1/(2*dt)*Cg);
updateStiffness=false;
Ug=gpuArray(U);

for i=2:1000
    Bstar=A1*Xi+A2*Xi1+Ug(i)*Forces;
    Xnew=A0\Bstar;
    X(:,i+1)=Xnew;
    Xi1=Xi;
    Xi=Xnew;
end

X=gather(X);
figure
yend=X(end-3,:);
plot(yend);

figure
thisframeNodes=Nodes;
animate=true;

if animate==true
   videoout=VideoWriter('HockeyStickNoPreload');
   videoout.FrameRate=30;
   open(videoout);
    for i=1:20
        j=i*50
        nodedisp=X(:,j);
        [thisframeNodes,thisframeCS_Origin]=newNodes(Nodes,nodedisp,Elements,DOF_Legend);
        
        
        PlotInitial(thisframeNodes,Elements,thisframeCS_Origin,false);
        drawnow;
     mframe=getframe;
     writeVideo(videoout,mframe);
     
    end
    close(videoout);
    
    
    %videoout.FileFormat='mp4';

%     writeVideo(videoout,movie(mframe,1,120))
    
    close(videoout);
else
end

EigenValues=diag(D);
[sortedEV,EVindex]=sort(abs(EigenValues));

% for i=1:2:9
% Vector1=V(n_dof+1:end,EVindex(i));
% Vector1=Vector1./max(Vector1).*0.5;
% 
%     [thisframeNodes,thisframeCS_Origin]=newNodes(Nodes,Vector1,Elements,DOF_Legend);
% 
%     PlotInitial(thisframeNodes,Elements,thisframeCS_Origin,false);
%     hold on
% end

% for i=1:length(EVindex)
%     ScaledOutput(i)=V(end-3,i)/max(V(n_dof+1:end,i));
% end
% ScaledOutput=ScaledOutput';
% [sortedOutput,outputindex]=sort(ScaledOutput);
% EVout=EigenValues(outputindex(end-5:end))./(2*pi);
