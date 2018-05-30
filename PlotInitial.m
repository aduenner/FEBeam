function [El_Origins,BodyCoords]=PlotInitial(Nodes,Elements,CS_Origin,CSdisplay,axisrange)

%CS_Origin is hte origin coordinate system of the elements which comes from
%getElementOrientation

CS0=[1 0 0 0;...
     0 1 0 0;...
     0 0 1 0;...
     0 0 0 1];
 
[nelements,~]=size(Elements);


for i=1:nelements
    %Determine HTM from origin of each element to global CS

   
   %Plot Original Position of Bodies                 
   nodefrom = Elements(i,2);
   nodeto   = Elements(i,3);
 
   bodyx(i)   = Nodes(nodefrom,2);
   bodyx(i+1)   = Nodes(nodeto,2);
      
   bodyy(i)   = Nodes(nodefrom,3);
   bodyy(i+1)   = Nodes(nodeto,3);
      
   bodyz(i)   = Nodes(nodefrom,4);
   bodyz(i+1)   = Nodes(nodeto,4);
   
   %Plot body
   

end
    ymin=min(bodyy)-.1;
   bodyplot=plot3(bodyx,bodyy,bodyz,'color','k','LineWidth',2);
if CSdisplay==true
    drawCS(CS_Origin);
else
end
BodyXYZ=[bodyx,bodyy,bodyz];
El_Origins = CS_Origin;
BodyCoords = BodyXYZ;
hold off

if nargin>4
    axis(axisrange)
else
    
axis([-0.05 0.5 ymin 0.1 -0.2 0.2])
end
daspect([1 1 1]);
view([135 -41])
set(gcf,'color','w');
% Set up figure properties:
% Enlarge figure to full screen.
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Get rid of tool bar and pulldown menus that are along top of figure.
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
