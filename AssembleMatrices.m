function [DOF_Legend,K,M,C] = AssembleMatrices(DOFs,Elements,ElementProperties,CS_Origin,alpha,beta)

tol=1e-6;
 
[nelements,~]=size(Elements);

%Initialize K as 12x12xn Vector and generate stiffness using BeamStiffness
K=zeros(12,12,nelements);
Mlocal=zeros(12,12,nelements);

for i=1:nelements
    ElementOrigins=CS_Origin(:,:,i);
    thisElementProperties=ElementProperties(i,:);
    K(:,:,i)=BeamStiffness(thisElementProperties,ElementOrigins);
    Mlocal(:,:,i)=BeamMass(thisElementProperties,ElementOrigins);
end


%Identify Unknown DOFs

%Global_DOF_Counter Count of unique dofs
%
%Global_DOF         Unique DOF, associated node and direction WRT node
%                   Col1=DOF, Col2=Node, Col3=direction 
%                   Col3 (1=x,2=y,3=z,4=rx,5=ry,6=rz)  
%
%RemainingDOF       Free DOFs in the local element k matrix
%                   1-6 are direction for node i and 7-12 for node j


Global_DOF_Counter=0;
for i=1:nelements
    
    %origin node of current element
    nodei=Elements(i,2); 
    
    %attachment node of current element
    nodej=Elements(i,3); 
    
    %Determine free DOF in local CS
    DOF_Mask=[DOFs(nodei,2:7),DOFs(nodej,2:7)]; %0=fixed 1=dof
    RemainingDOF=(1:1:12).*DOF_Mask;
    RemainingDOF(RemainingDOF==0)=[];
    
    %RemainingDOF has up to 12DOF. 1-6-->node1 and 7-12-->node2
    for j=1:length(RemainingDOF)
            this_dof=RemainingDOF(j);
        if this_dof>6
            this_dof=this_dof-6;
            participating_node=nodej;
        else
            participating_node=nodei;
        end
        
        %Potentially add new DOF to global unknown DOF
        %this_dof is the local dof 1=x 2=y 3=z 4=rx 5=ry 6=rz
        %of participating_node element i
        
        Append_dof_legend_vals=[participating_node,this_dof];
        
        %Add first row to Global_DOF if it is empty
        if Global_DOF_Counter==0
        Global_DOF=[1,Append_dof_legend_vals];
        Global_DOF_Counter=1;
        end
        
        %For rows after the first row, check for duplicates of node and
        %node dof direction. Keep track of rows where the duplicates occur
        %so that we can associate nodal DOFs in element stiffness matrices
        %with the global list of DOFs

        if Global_DOF_Counter>=1
            
            %Loop through all Global_DOF rows
            for i1=1:Global_DOF_Counter
                doftocheck=Global_DOF(i1,2:3);
                duplicatecheck=abs(doftocheck-Append_dof_legend_vals);
                
                %Check row for duplicate values
                if sum(duplicatecheck)<tol
                    duplicaterow=i1;
                    
                    %Store row for later use in matching local/global DOF
                    Alegend(j)=duplicaterow;
                    
                    %Exit loop if there is a duplicate, move to next
                    %RemainingDOF value (j)
                    break
                
                else
                    %If no duplicates found after reaching the end of the 
                    %unique DOF list, we do the following
                    if i1==Global_DOF_Counter
                        Global_DOF_Counter=Global_DOF_Counter+1;
                        Global_DOF=[Global_DOF;[Global_DOF_Counter,Append_dof_legend_vals]];
                        Alegend(j)=Global_DOF_Counter;
                    else
                    end
                end 
            end
            
%         else
%             duplicaterow=0;
%             Global_DOF_Counter=Global_DOF_Counter+1;
%             Global_DOF=[Global_DOF;[Global_DOF_Counter,Append_dof_legend_vals]];
        end

    end
    %Record the row of unique global unknown DOFS for association with the 
    %local stiffness matrix of each element
    Local_k_legend{i}=Alegend;
    Local_m_legend{i}=Alegend;
    
    %Reduce stiffness matrix based on unique DOF
    K_orig=K(:,:,i);
    K_reduced{i}=K_orig(RemainingDOF,RemainingDOF);
    
    M_orig=Mlocal(:,:,i);
    M_reduced{i}=M_orig(RemainingDOF,RemainingDOF);
    
end


GlobalStiffnessDOF=Global_DOF_Counter;

%GlobalStiffness (reduced stiffness matrix) is of size nxn where n=
%                number of unique DOFs. In the next section we assemble it

GlobalStiffness = zeros(GlobalStiffnessDOF);
     GlobalMass = zeros(GlobalStiffnessDOF);

%Col1=Global_DOF_counter value, Col2=Node, Col3=local dof (x=1,y=2,etc)
%Keep track of unknown DOFs added to reduced stiffness matrices

for globalrow=1:GlobalStiffnessDOF
    
    %for each row in GlobalStiffness (unique dof)

    %Search local reduced stiffness matrices of each element
    for element_iter=1:nelements
        %pull local stiffness matrix for current element
        thisK=K_reduced{element_iter};
        thisM=M_reduced{element_iter};
        
        thisk_legend=Local_k_legend{element_iter};
        thism_legend=Local_m_legend{element_iter};
        
        %determine #dof of current element stiffness matrix
        [thisK_DOF,~]=size(thisK);
        [thisM_DOF,~]=size(thisM);
            
         %IF 
            %local dof represented by row of component of local stiffness 
            %matrix matches dof represented by current row in global matrix
            
        %THEN
            %Add local stiffness component from each column in local
            %element stiffness matrix to the matching column in the current
            %row of the global stiffness matrix
            
            for DOF=1:thisK_DOF
                if thisk_legend(DOF)==globalrow
                    for iter1=1:thisK_DOF
                        Kij=thisK(DOF,iter1);
                        Mij=thisM(DOF,iter1);
                        globalstiffnesscol=thisk_legend(iter1);
                             globalmasscol=thism_legend(iter1);
                             
                        GlobalStiffness(globalrow,globalstiffnesscol)=...
                        GlobalStiffness(globalrow,globalstiffnesscol)+Kij;
                    
                        GlobalMass(globalrow,globalstiffnesscol)=...
                            GlobalMass(globalrow,globalstiffnesscol)+Mij;
                    
                    end
                end
            end
    end
end

K=GlobalStiffness;
M=GlobalMass;
DOF_Legend=Global_DOF;

C=alpha*M+beta*K;