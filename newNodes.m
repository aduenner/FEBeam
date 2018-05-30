function [postNodes,postElement_CS]=newNodes(Nodes,NodalDisplacements,Elements,DOF_Legend)

[n_dof,~]=size(DOF_Legend);
[n_el,~]=size(Elements);

Nodes_post=Nodes;

for i=1:n_dof
    active_node=DOF_Legend(i,2);
    dof_index=DOF_Legend(i,3);
    if dof_index<7
    Nodes_post(active_node,dof_index+1)=Nodes(active_node,dof_index+1)+...
                                      NodalDisplacements(i);
    end
end
postNodes=Nodes_post;
postElement_CS=getElementOrientation(Nodes_post,Elements);
