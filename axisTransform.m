function Tmatrix = axisTransform(angles,displacements)

rx=angles(1);
ry=angles(2);
rz=angles(3);
x=displacements(1);
y=displacements(2);
z=displacements(3);

%rx=Rotation about this CS X axis to align y axes
%ry=Rotation about resulting Y axis to align axes
%rz=Rotation about resulting z axis to align axes
%tx=Translation about new x axis
%ty=Translation about new y axis
%ty=Translation about new z axis

%Note: Transformation occurs in order CS1-->rx-->ry-->rz-->txtytz-->CS2

rotx=[1     0           0          0;...
      0     cos(rx)     -sin(rx)   0;...
      0     sin(rx)     cos(rx)    0;...
      0     0           0          1];
  
  
roty=[cos(ry)   0   sin(ry)     0;...
      0         1   0           0;...
      -sin(ry)  0   cos(ry)     0;...
      0         0   0           1];
  
  rotz=[cos(rz)  -sin(rz)   0   0;...
        sin(rz)   cos(rz)   0   0;...
        0         0         1   0;...
        0         0         0   1];   
    
    trans=[1 0 0 x;...
           0 1 0 y;...
           0 0 1 z;...
           0 0 0 1];
    

    
    Tmatrix=rotx*roty*rotz*trans;
