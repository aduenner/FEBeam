function M = BeamMass(ElementProperties,Orientation)

%Input ElementProperties,Orientation
    %ElementProperties -- Vector containing [E,rho,nu,A,Iy,Iz,J,L,G]
    %Orientation       -- 4x4 HTM from (CS_Origin of member) to Global_CS; 

E   = ElementProperties(1);     %Elastic Modulus
rho = ElementProperties(2);     %Density
nu  = ElementProperties(3);     %Poisson's Ratio
A   = ElementProperties(4);     %Cross-Sectional Area
Iy  = ElementProperties(5);     %Moment of Inertia about Y
Iz  = ElementProperties(6);     %Moment of Inertia about Z
J   = ElementProperties(7);     %Polar Moment of Inertia 
L   = ElementProperties(8);     %Beam Length

%Get rotation angles from the HTM
RotationAngles=Orientation(1:3,1:3);

%Generate block diagonal 12x12 matrix w/ 'rotationblock' on diagonals
RotationTransform=blkdiag(RotationAngles,RotationAngles,RotationAngles,RotationAngles);

%Assembly 4 6x6 Matrices 
%K=[K11 K12
%  [K21 K22] 
%where Kij is compliance matrix from node i to node j
I0=J;


Mlocal=zeros(12,12);
Mlocal(12,1:11) = [0  -13*L 0     0        0      -3*L^2 0  -22*L 0    0 0];
Mlocal(11,1:10) = [0  0     13*L  0        -3*L^2 0      0  0     22*L 0];
Mlocal(10,1:9)  = [0  0     0     70*I0/A  0      0      0  0     0];
Mlocal(9,1:8)   = [0  0     54    0        -13*L  0      0  0];
Mlocal(8,1:7)   = [0  54    0     0        0      13*L   0];
Mlocal(7,1:6)   = [70 0     0     0        0      0];
Mlocal(6,1:5)   = [0  22*L  0     0        0];
Mlocal(5,1:4)   = [0  0     -22*L 0];
Mlocal(4,1:3)   = [0  0     0];
Mlocal(3,1:2)  = [0  0];
Mlocal(2,1)    =  0;

Mdiag=diag([140 156 156 140*I0/A 4*L^2 4*L^2 140 156 156 140*I0/A 4*L^2 4*L^2]);
Mlocal=Mlocal+Mdiag;
Mlocalsymm=Mlocal+tril(Mlocal,-1).';

MLocalScaled=Mlocalsymm.*rho*A*L/420;

 
 %Output K in Global Coordinates
M=RotationTransform*MLocalScaled*RotationTransform';
