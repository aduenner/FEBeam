function [sortedEV,EVindex]=  HockeyEigen(A,Nodes,Elements,DOF_Legend,index)

tic
[V,D]=eig(A);
EigenValues=diag(D);
[sortedEV,EVindex]=sort(abs(EigenValues));
eigentime=toc
[n_dof,~]=size(A);
Vector1=V((end/2+1):end,EVindex(index));
Vector1=real(Vector1);
scalefactor=1/max(abs(Vector1));
Vector1scaled=Vector1*scalefactor;
maxm=1;
scale=[linspace(0,maxm,15) linspace(maxm,0,15) linspace(0,-maxm,15) linspace(-maxm,0,15)];
for f=1:60
    mfactor=scale(f)/max(abs(Vector1scaled(3:3:end)))*0.2;
    EVscaled=Vector1scaled.*mfactor;
    
    [thisframeNodes,thisframeCS_Origin]=newNodes(Nodes,EVscaled,Elements,DOF_Legend);

    PlotInitial(thisframeNodes,Elements,thisframeCS_Origin,false);
    %.2f
    titlestr=['Eigenfrequency = ',num2str(sortedEV(index)./(2*pi)),'hz'];
    title(titlestr);
    thisframeNodes(end-3);
    drawnow
end
% 
% for i=1:length(EVindex)
%     ScaledOutput(i)=V(end-3,i)/max(V(n_dof+1:end,i));
% end
% ScaledOutput=ScaledOutput';
% [sortedOutput,outputindex]=sort(ScaledOutput);
% EVout=EigenValues(outputindex(end-5:end))./(2*pi);