%drawCS.mat

%Function - Plot coordinate system axes
%Input - 4x4xn array of homogeneous coordinate systems WRT global CS
%Output - graph w/ 3 orthogonal unit vectors aligned with local cs u,v,w
%         located at origin of each local coordinate system


%Andrew Duenner
%November 21, 2016


function drawCS(CSmatrices)
%Scale quiver plots appropriately
scale=0.05;
originoffset=scale/2;
%Get number of coordinate systems
[~,~,numCS]=size(CSmatrices);
texton=false;
hold on
for i=1:numCS
    pose=CSmatrices(:,:,i);
    CSnumber=i;
    
    origin=pose(1:3,4)';
    
    udir=pose(1:3,1)';
    vdir=pose(1:3,2)';
    wdir=pose(1:3,3)';
    
    udirloc=origin+scale*udir;
    vdirloc=origin+scale*vdir;
    wdirloc=origin+scale*wdir;
    
    origintext=[num2str(i-1)];
    %origintext=origintext{:};
    udirtext={'X'};
    vdirtext={'Y'};
    wdirtext={'Z'};
    
    nulldir=[0,0,0];
    
    uquiver=num2cell([origin(:);udir(:)]');
    vquiver=num2cell([origin(:);vdir(:)]');
    wquiver=num2cell([origin(:);wdir(:)]');
    quiver3(uquiver{:},scale,'r')
    quiver3(vquiver{:},scale,'g')
    quiver3(wquiver{:},scale,'b')
    
    
    origintextloc=num2cell((origin-[originoffset originoffset originoffset]));
    utextloc=num2cell(udirloc(:));
    vtextloc=num2cell(vdirloc(:));
    wtextloc=num2cell(wdirloc(:));
    
    if texton==true
        text(origintextloc{:},origintext);
        text(utextloc{:},udirtext);
        text(vtextloc{:},vdirtext);
        text(wtextloc{:},wdirtext);
    else
    end
end
