 function K = BeamStiffness(ElementProperties,Orientation)
 
%Input ElementProperties,Orientation
    %ElementProperties -- Vector containing [E,nu,A,Iy,Iz,J,L,G]
    %Orientation       -- 4x4 HTM from (CS_Origin of member) to Global_CS; 

E  = ElementProperties(1);     %Elastic Modulus
nu = ElementProperties(3);     %Poisson's Ratio
A  = ElementProperties(4);     %Cross-Sectional Area
Iy = ElementProperties(5);     %Moment of Inertia about Y
Iz = ElementProperties(6);     %Moment of Inertia about Z
J  = ElementProperties(7);     %Polar Moment of Inertia 
L  = ElementProperties(8);     %Beam Length

%Get rotation angles from the HTM
RotationAngles=Orientation(1:3,1:3);

%Generate block diagonal 12x12 matrix w/ 'rotationblock' on diagonals
RotationTransform=blkdiag(RotationAngles,RotationAngles,RotationAngles,RotationAngles);

%Assembly 4 6x6 Matrices 
%K=[K11 K12
%  [K21 K22] 
%where Kij is compliance matrix from node i to node j

K11=diag([A/L, 12*Iz/L^3, 12*Iy/L^3,J/(2*(1+nu)*L),...
          4*Iy/L,4*Iz/L]);
K11(6,2)=6*Iz/L^2;
K11(5,3)=-6*Iy/L^2;
K11(3,5)=-6*Iy/L^2;
K11(2,6)=6*Iz/L^2;

K12=diag([-A/L, -12*Iz/L^3, -12*Iy/L^3,-J/(2*(1+nu)*L),...
          2*Iy/L,2*Iz/L]);
K12(6,2)=-6*Iz/L^2;
K12(5,3)=6*Iy/L^2;
K12(3,5)=-6*Iy/L^2;
K12(2,6)=6*Iz/L^2;

K21=diag([-A/L, -12*Iz/L^3, -12*Iy/L^3,-J/(2*(1+nu)*L),...
          2*Iy/L,2*Iz/L]);
K21(6,2)=6*Iz/L^2;
K21(5,3)=-6*Iy/L^2;
K21(3,5)=6*Iy/L^2;
K21(2,6)=-6*Iz/L^2;

K22=diag([A/L, 12*Iz/L^3, 12*Iy/L^3,J/(2*(1+nu)*L),...
          4*Iy/L,4*Iz/L]);
K22(6,2)=-6*Iz/L^2;
K22(5,3)=6*Iy/L^2;
K22(3,5)=6*Iy/L^2;
K22(2,6)=-6*Iz/L^2;

Klocal=E*[K11,K12;...
     K21,K22];
 
 %Output K in Global Coordinates
K=RotationTransform*Klocal*RotationTransform';
