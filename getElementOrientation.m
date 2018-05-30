
function orientationOut = getElementOrientation(Nodes,Elements)

%Get number of elements and nodes
[numel,~]=size(Elements);
[numnodes,~]=size(Nodes);
orientationOut=zeros(4,4,numel);
for i=1:numel
    %Origin node is in 2nd col of element aray
    %Destination node is in 3rd col of element array
    node1=Elements(i,2);
    node2=Elements(i,3);
    
    element_angle=Nodes(node1,5:7);
    
    %Get pose of element's origin node 
    element_origin_pose=axisTransform([0,0,0],[0,0,0])*...
                        axisTransform(element_angle,[0,0,0]);
    
    element_y=element_origin_pose(1:3,2)';
    
    %XYZ coordinates located in cols 2-4 of node array
    node1xyz=Nodes(node1,2:4);
    node2xyz=Nodes(node2,2:4);
    
    %First direction cosine in direction from origin to destination node
    i1=(node2xyz-node1xyz)/norm(node2xyz-node1xyz);
    
    %Third point is in direction of yhat in origin node
    p_xyz=i1+element_y;
    
    %Vectors from n1 to n2
    v21=node2xyz-node1xyz;
    v31=p_xyz-node1xyz;
    
    l=norm(v21);
    xhat=v21/l;
    zhat=cross(v21,v31)/norm(cross(v21,v31));
    yhat=cross(zhat,xhat)/norm(cross(zhat,xhat));

i2=cross([0 0 1],i1);
i3=cross(i1,i2);


orientationOut(:,:,i)=[[[i1;i2;i3]';[0 0 1]] [node1xyz';0]];
end
