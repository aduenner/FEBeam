function [Nodes,Nodal_DOFS,Elements,El_Properties]=HockeyInitialNodes(ShaftNodes,BladeNodes)



%Geometry
ShaftLength=0.7;
ShaftH=1e-2;
ShaftW=3e-2;

BladeLength=.3;
BladeWidth=1e-2;
BladeTH=3e-2;
BladeAngle=-pi/2+10*pi/180;

%%SECTION PROPERTIES
E=7e10;
rho=7000;
Nu=0.33;

Ashaft=ShaftH*ShaftW;
Iyshaft=ShaftW*ShaftH^3/12;
Izshaft=ShaftH*ShaftW^3/12;
Jshaft=ShaftH*ShaftW*(ShaftH^2+ShaftW^2)/12;

Ablade=BladeWidth*BladeTH;
Iyblade=BladeWidth*BladeTH^3/12;
Izblade=BladeTH*BladeWidth^3/12;
Jblade=BladeWidth*BladeTH*(BladeWidth^2+BladeTH^2)/12;

dL_shaft=ShaftLength/(ShaftNodes);
dL_blade=BladeLength/(BladeNodes);

El_Properties_Shaft=[E rho Nu Ashaft Iyshaft Izshaft Jshaft dL_shaft];
El_Properties_Blade=[E rho Nu Ablade Iyblade Izblade Jblade dL_blade];


%Initialize Arrays

% id, x  y  z  rx  ry      rz
Nodes=[1, 0, 0, 0, 0,  -pi/2,  0];
% id  x  y  z  rx  ry  rz
Nodal_DOFS=[1, 0, 0, 0, 0,   0,     0];

Elements=[1 1 2];

El_Properties=El_Properties_Shaft;

%%Generate Stick

nextnode=Nodes(end,1)+1;
for i=nextnode:ShaftNodes
    Nodes(i,:)=[i 0 -dL_shaft*(i-1) 0 0 -pi/2 0];
    Nodal_DOFS(i,:)=[i 1 1 1 1 1 1];
    if i>2
        Elements(i-1,:)=[i-1 Nodes(i-1,1) Nodes(i,1)];
        El_Properties(i-1,:)=El_Properties_Shaft;
    else
    end
end

%%Generate Transition

% Nodes(end,6)=-BladeAngle;
nextnode=Nodes(end,1)+1;
for i=nextnode:nextnode+BladeNodes
    EndPos=Nodes(i-1,2:3);
    
    if i==nextnode
        yTransition=-dL_shaft;
        xTransition=0;
        El_Properties(end+1,:)=El_Properties_Shaft;
    else
        yTransition=dL_blade*sin(BladeAngle);
        xTransition=dL_blade*cos(BladeAngle);
        El_Properties(end+1,:)=El_Properties_Blade;
    end
    
    FaceNodeEnd=EndPos+[xTransition,yTransition];
    
    Nodes(i,:)=[i,FaceNodeEnd(:)',0,0,-BladeAngle,0];
    Nodal_DOFS(i,:)=[i 1 1 1 1 1 1];
    Elements(end+1,:)=[Elements(end,1)+1,i-1,i];
    
end

