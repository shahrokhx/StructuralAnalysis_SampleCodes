%%                       P L A N A R     F R A M E                  
%__________________________________________________________________________
% 
%                          Structural Analysis
%                            Dr. Rafi Muhanna
%                   Georgia Institute of Technology
%                              Spring 2019
%__________________________________________________________________________
%  Developed by 
%  SHAHROKH SHAHI  (www.sshahi.com)
%  (based on initial version by ...)
%
%  Usage: 
%  The only inputs that you need to provide is:
%  (1) input file name (LINE 28) 
%  (2) output file name (LINE 29)
%  (3) a (0/1) switch value for displaying stiffness matrices (LINE 32) 
% 
%% Initialization
clc
clear 
close all

format short g
format compact

%% Inputs
inpFileName = 'input.txt';
outFileName = 'output.txt';

% Do you want to see the stiffness matrices?  (1 = Yes,  0 = No)
showStiff   = 1; 

%% Body of the Procedure

% setting up the output file
outFile = fopen(outFileName,'w');

% reading input file
Model = inpFileReader(inpFileName);

% displaying summary of input file on Command Window and output file
printSummary(Model)
printSummary(Model, outFile)

% graphical view of the geometry (on figure #1)
PlanePlot(Model,1)

% analysing the model
Solution = frameAnalysis(Model);

% displaying siffness matrices (if needed)
if showStiff
    printStiff(Solution)
    printStiff(Solution, outFile)
end

% displaying a summary of results on Command Window and output file
printResult(Solution)
printResult(Solution, outFile)

% displaying deformed shape
figureNumber = 2;
PlanePlotDeform(Model, Solution, figureNumber)

% closing the output file
fclose(outFile);


%% Functions
%--------------------------------------------------------------------------
function Solution = frameAnalysis(Model)

% Model = ModelRevise(Model);

Solution.info = ['Solution obtained on ', date];

FirstNode  = Model.geometry.elements(:,1);
SecondNode = Model.geometry.elements(:,2);
ifix       = Model.bound.ifix;

% processing element load information
ndof = Model.nNode * Model.nDof;
K    = zeros(ndof,ndof); %initialize the global stiffness matrix of size ndof
Feq  = zeros(ndof,1); %nodal force vector due to element loads
ke_local (:,:,Model.nElem) = zeros(6);
ke_global(:,:,Model.nElem) = zeros(6);
T6(:,:,Model.nElem)        = zeros(6);
dofList                    = zeros(Model.nElem,6); % to store the assembly indices

for e = 1 : Model.nElem
    c  = Model.geometry.Cx(e);
    s  = Model.geometry.Cy(e);
    le = Model.geometry.L(e);

    ae = Model.properties.A(e) * Model.properties.E(e); %axial rigidity 
    ei = Model.properties.E(e) * Model.properties.Iz(e); %flexural rigidity
    %Degrees of freedom for element e with nodes m and n are 3m-2, 3m-1 ,3m, 3n-2, 3n-1 ,3n
    %These degrees of freedom are stored in the vector {dof} as shown below
    dof = [3*FirstNode(e)-2, 3*FirstNode(e)-1, 3*FirstNode(e), ...
           3*SecondNode(e)-2, 3*SecondNode(e)-1, 3*SecondNode(e)];
       
    % rotation transformation matrix for element e
    T=[c    s    0    0     0     0;
      -s    c    0    0     0     0;
       0    0    1    0     0     0;
       0    0    0    c     s     0;
       0    0    0   -s     c     0;
       0    0    0    0     0     1;];
   
    % element stiffness matrix w.r.t local axes
    estif = [ ae/le     0             0          -ae/le     0            0; 
               0       12*ei/le^3    6*ei/le^2       0    -12*ei/le^3   6*ei/le^2;
               0       6*ei/le^2     4*ei/le         0    -6*ei/le^2    2*ei/le ;
            -ae/le      0             0           ae/le      0           0; 
               0      -12*ei/le^3   -6*ei/le^2       0     12*ei/le^3  -6*ei/le^2;
               0       6*ei/le^2     2*ei/le         0    -6*ei/le^2    4*ei/le ];
               
    if Model.geometry.hingIndex(SecondNode(e)) == 1   % if the right node is a hinge
        % Modifications due to the hing on the right hand
        estif( 6   , :   ) = 0;
        estif( :   , 6   ) = 0;
        estif([2,5],[2,5]) = (1/4) * estif([2,5],[2,5]);
        estif([2,5], 3   ) = (1/2) * estif([2,5], 3   );
        estif( 3   ,[2,5]) = (1/2) * estif( 3   ,[2,5]);
        estif( 3   , 3   ) = (3/4) * estif( 3   , 3   );
    end
     
     
     ke=(T'*estif*T) ;  %element stiffness matrix with reference to global axes
     
     ke_local(:,:,e) = estif;
     T6(:,:,e)= T;
     ke_global(:,:,e) = ke;
     dofList(e,:) = dof;
     
     K(dof,dof)=K(dof,dof)+ke;  %assembly of global stiffness matrix K
     %{Feq} is vector of equivalent nodal forces due to element loads
     Feq(dof) = Feq(dof) + T'*Model.loading.eforce(:,e);     
end

Fc = Model.loading.F+Feq; %combined force vector {Fc} is obtained by adding vector of nodal loads {F} to {Feq}

Solution.data.ke_local = ke_local;
Solution.data.ke_global = ke_global;
Solution.data.T6 = T6;
Solution.data.F_global = Fc;
Solution.data.dofList = dofList;
Solution.data.K = K; % assembled global stiff matrix
%-------------------------------------------------------------------------
%applying boundary conditions (specified by restraints along X and Y axes)
%on global stiffness matrix and global force vector
K(ifix,:)=0; %Elements in rows with d.o.f corresponding to ifix are made 0                 
K(:,ifix)=0; %Elements in columns with d.o.f corresponding to ifix are made 0                    
K(ifix,ifix)=eye(length(ifix));%Diagonal elements with d.o.f corresponding to ifix are made 1     
Fc(ifix) = 0; %corresponding columns of Global force vector are also made zero
%---------------------------------------------------------------------
%obtaining vector of nodal displacements along global axes
u=inv(K)*Fc;

Solution.data.u = u;


df = setdiff(1:length(Fc),ifix);
Solution.data.K_bc = K(df,df);



%---------------------------------------------------------------------
%P O S T P R O C E S S I N G
%---------------------------------------------------------------------
%Computation of internal forces and moments in elements of Plane Frame
for e = 1 : Model.nElem
    %degrees of freedom corresponding to element e
    dof=[3*FirstNode(e)-2 3*FirstNode(e)-1 3*FirstNode(e) 3*SecondNode(e)-2 3*SecondNode(e)-1 3*SecondNode(e)];

    % Vector of internal forces of element w.r.t. local axes   
    % lforce(:,e) = estif*T*u(dof)-Model.loading.eforce(:,e);   
    lforce(:,e) = ke_local(:,:,e) * T6(:,:,e) * u(dof)-Model.loading.eforce(:,e);
end

Solution.data.lforce = lforce;
Solution.nodalDisplacements = reshape(u,Model.nDof,Model.nNode)';


% Calculating the Reactions
reactionForces = [];
for i = 1 : Model.nNode
    if(Model.bound.resnode(i)==0) %restrained nodes (at support) are identified 
       Rx = 0.0; %components of support reactions along global X and Y directions are initialized
       Ry = 0.0; 
       Mz = 0.0;%support moment is initialized
       for e = 1 : Model.nElem %elements with "node i" as one of the nodes are identified
          %components of force in element e along global X and Y
          %directions are extracted using direction cosines and are
          %then its contribution to Rx and Ry is computed
           if(i == FirstNode(e))  %"node i" is the first node of element e
              Rx = Rx + Model.geometry.Cx(e)*lforce(1,e) - Model.geometry.Cy(e)*lforce(2,e);
              Ry = Ry + Model.geometry.Cy(e)*lforce(1,e) + Model.geometry.Cx(e)*lforce(2,e);
              Mz = Mz + lforce(3,e);
           end
           if(i == SecondNode(e)) %"node i" is the second node of element e
              Rx = Rx + Model.geometry.Cx(e)*lforce(4,e) - Model.geometry.Cy(e)*lforce(5,e);
              Ry = Ry + Model.geometry.Cy(e)*lforce(4,e) + Model.geometry.Cx(e)*lforce(5,e);
              Mz = Mz + lforce(6,e);
           end
       end
       %support reactions are set zero for any unrestrained degree of freedom
       %associated with a support
       if(Model.bound.resx(i) == 1)  
           Rx = 0.0;
       end
       if(Model.bound.resy(i) == 1)
           Ry = 0.0;
       end
       if(Model.bound.resz(i) == 1)
           Mz = 0.0;
       end
       reactionForces = [reactionForces; i,Rx,Ry,Mz];
    end
end
Solution.reactionForces = reactionForces;






%--------------------------------------------------------------------------
% this part will be moved in the future
if isfield(Model,'analysis')
    if Model.analysis.showResults
        printResult(Solution,1)
    end
    if Model.analysis.showDetails
        % print stiffness matrices
        printStiff(Solution,1)
    end
    if Model.analysis.saveToFile
        outFile = fopen(Model.analysis.outputFileName, 'w');
        printSummary(Model,outFile);
        printResult(Solution,outFile);

        fclose(outFile);
    end
end
%--------------------------------------------------------------------------
end
%--------------------------------------------------------------------------
function Model = inpFileReader(inpFileName)
inp   = fopen(inpFileName,'r');
%--------------------------------------------------------------------------
%       R E A D I N G     G E N E R A L      I N F O R M A T I O N
%--------------------------------------------------------------------------
mat   = fscanf(inp,'%d',2); % this line reads number of nodes and number of 
                            % elements in the truss
nNode = mat(1);             % number of nodes 
nElem = mat(2);             % number of elements

nDof  = 3*nNode;            % total number of DOFs for the Plane Frame

%--------------------------------------------------------------------------
%         R E A D I N G      N O D A L      I N F O R M A T I O N
%--------------------------------------------------------------------------

%Read nodal information from the input file node-by-node
mat   = fscanf(inp,'%d %f %f %d %d %d  %e %e  %e  %d',[10,nNode]);

mat   = mat';     % transposed for future use
Xval  = mat(:,2); % x coordinates of nodes (column 2 of input as above)
Yval  = mat(:,3); % y coordinates of nodes (column 3 of input as above)

% restraints of a node
resx  = mat(:,4); % restraint along x-direction from column 4
resy  = mat(:,5); % restraint along y-direction from column 5 
resz  = mat(:,6); % restraint about z-direction from column 6 

%--------------------------------------------------------------------------
%  Generation of  nodal force vector (F) (of size ndof) along global axes
%--------------------------------------------------------------------------
% degrees of freedom corresponding to i'th node are 3i-2, 3i-1 and 3i
% fx , fy  and Mz are components of nodal forces along global X and Y 
% directions and Moment about z-direction respectively
% At node i, fx acts along degree of freedom  '3i-2' 
%            fy acts along degree of freedom  '3i-1' 
%            Mz acts about degree of freedom  '3i' 
fx   = mat(:,7);      % reading fx from column 7
fy   = mat(:,8);      % reading fy from column 8
Mz   = mat(:,9);      % reading Mz from column 9

hingeIndex = mat(:,10);  % hinge = 1    regular node = 0
F    = zeros(nDof,1); % initializing global nodal force vector
for i = 1 : nNode
  F(3*i-2) = fx(i);  
  F(3*i-1) = fy(i); 
  F(3*i)   = Mz(i); 
end

%--------------------------------------------------------------------------
%      B O U N D A R Y    C O N D I T I O N    I N F O R M A T I O N
%--------------------------------------------------------------------------
% bn is a counter for the total number of restraints
bn = 0; % initialized to zero
%--------------------------------------------------------------------------
% notation used to indicate restraint associated with a given dof
% 0 indicates RESTRAINED   degree of freedom
% 1 indicates UNRESTRAINED degree of freedom
%--------------------------------------------------------------------------
% A new vector "{resnode}" of size nNode is initialized with ones 
% {resnode} is used to indicate the status of restraint for nodes:
%  - Restrained   nodes (at supports): 0  
%                       (with atleast one RESTRAINED degree of freedom) 
%  - Unrestrained nodes : 1resnode = ones(nNode,1);

resnode = ones(nNode,1); % nodes are initialized as UNRESTRAINED by default

for i = 1 : nNode
   if(resx(i) == 0) % if this dof along x-direction is restrained
        bn = bn+1;
        ifix(bn) = 3*i-2;  % add to the list of restrained dofs
        resnode(i) = 0; 
   end
   if(resy(i) == 0) % if this dof along x-direction is restrained
        bn = bn+1;
        ifix(bn) = 3*i-1;  % add to the list of restrained dofs
        resnode(i) = 0;
   end
   if(resz(i) == 0) % if this dof ABOUT z-direction is restrained
        bn = bn+1;
        ifix(bn) = 3*i;    % add to the list of restrained dofs
        resnode(i) = 0;
   end
end
% NOTE: bn is incremented by 3  for a fixed support 
%       bn is incremented by 2  for a pin support 
%       bn is incremented by 1  for a roller support in the above for-loop
% at the end of the above loop, 
%       bn   is total number of restrained degrees of freedom
%       ifix is list of restrained degrees of freedom for the Plane Frame

%--------------------------------------------------------------------------
%  R E A D I N G     E L E M E N T     I N F O R M A T I O N
%--------------------------------------------------------------------------
mat = fscanf(inp,'%d %d %d %f %e %e %f %f %f %f',[10,nElem])';
% mat = mat'
Cx      = zeros(nElem,1);    % direction cosines of elements Cx and Cy
Cy      = zeros(nElem,1);
L       = zeros(nElem,1);    % length of each element
eforce  = zeros(6,nElem);    % element force vector with six dofs 
                             % (3 at each node of element)
                             
FirstNode  = mat(:,2);
SecondNode = mat(:,3);
A       = mat(:,4); % cross sectional areas of elements from column 4
E       = mat(:,5); % Modulus of Elasticity of elements from column 6
Iz      = mat(:,6); % area moment of inertia of elements from column 5

% reading element load information
% uniformly distributed load of intensity w1 acts from a distance of d1 to d2 
% as measured from left end of element
wr = mat(:,7);  % distributed load w1  acting on elements from column 7
d1 = mat(:,8);  % left limit d1 of distributed load w1 from columns 8
wl = mat(:,9);  % distributed load w2  acting on elements from column 9
d2 = mat(:,10); % right limit d2 of distributed load w1 from column 10
w1 = -wr;


%Computation of equivalent nodal loads owing to uniformly distributed load over some portion of the span
for e = 1 : nElem
    x1 = Xval(FirstNode(e));    
    y1 = Yval(FirstNode(e));  %x and y coordinates of first node
    x2 = Xval(SecondNode(e));   
    y2 = Yval(SecondNode(e)); %x and y coordinates of second node
    le = sqrt((x2-x1)^2+(y2-y1)^2); %length of element e
    %direction cosines of element e
    c  = (x2-x1)/le;
    s  = (y2-y1)/le;
    %Storing values of direction cosines c and s in vectors {Cx} and {Cy} for subsequent use 
    Cx(e) = c;
    Cy(e) = s;
    L (e) = le;
    %computation of equivalent nodal loads owing to uniformly distributed load over some portion of the span 
    c1 = d2(e) - d1(e);
    a  = d1(e) + 0.5*c1;
    b  = 0.5*c1 + le - d2(e);
    p1 = (-w1(e)*c1 / (12*le^2)) * (12*a*b^2 + c1^2*(le-3*b));
    p2 = ( w1(e)*c1 / (12*le^2)) * (12*a^2*b+c1^2*(le-3*a));
    p3 = ( w1(e)*c1 * b/le) - (p1+p2)/le;
    p4 = ( w1(e)*c1 * a/le) + (p1+p2)/le;
    eforce(:,e) = eforce(:,e) + [0.0 -p3 p1 0.0 -p4 p2]';
    %Negative of vector of fixed end forces in added to the element force vector
end
%--------------------------------------------------------------------------
%            F O R M I N G     M O D E L     S T R U C T U R E 
%--------------------------------------------------------------------------
Model.info         = 'Plane Beam/Frame Program';
Model.analysisType = 'BEAM';
Model.nDim         = 2;  % 2D (plane beam/frame)
Model.nNode        = nNode;
Model.nElem        = nElem;
Model.nDof         = 3;
Model.nElemNode    = 2;

Model.geometry.coordinates = [Xval, Yval];
Model.geometry.elements    = [FirstNode, SecondNode];
Model.geometry.Cx  = Cx;
Model.geometry.Cy  = Cy;
Model.geometry.L   = L;
Model.geometry.hingIndex = hingeIndex;

Model.properties.A = A;
Model.properties.E = E;
Model.properties.Iz= Iz;

Model.nBoundary    = bn;
Model.bound.resx = resx;
Model.bound.resy = resy;
Model.bound.resz = resz;
Model.bound.resnode = resnode;
Model.bound.ifix = ifix;

Model.loading.F = F;
Model.loading.eforce = eforce;
Model.loading.data = [d1,d2,wl,wr];
end
%--------------------------------------------------------------------------
function printSummary(Model,out)
if nargin < 2
    out = 1;
end
Xval = Model.geometry.coordinates(:,1);
Yval = Model.geometry.coordinates(:,2);
resx = Model.bound.resx;
resy = Model.bound.resy;
resz = Model.bound.resz;
F    = Model.loading.F;
eforce = Model.loading.eforce;
d1   = Model.loading.data(:,1);
d2   = Model.loading.data(:,2);
wl   = Model.loading.data(:,3);
wr   = Model.loading.data(:,4);

fprintf(out,'A N A L Y S I S   O F   A   P L A N E   F R A M E \n\n');
fprintf(out,'N O D A L   I N F O R M A T I O N \n');
fprintf(out,'________________________________________________________________________________________________\n');
fprintf(out,'Node   Coordinates      *Nodal          Nodal Loads    In-plane    Type of support\n');
fprintf(out,'                      Restraints                       moment\n');
fprintf(out,'------------------------------------------------------------------------------------------------\n');
fprintf(out,'       X       Y      x   y theta(z)    Fx      Fy       Mz\n');
fprintf(out,'------------------------------------------------------------------------------------------------\n');
for i = 1 : Model.nNode
     fprintf(out,'%2d   %4.1f    %4.1f     %d   %d   %d    %7.1f  %7.1f %7.1f',i,Xval(i),Yval(i),resx(i),resy(i),resz(i),F(3*i-2),F(3*i-1),F(3*i)); 
     if(resx(i)==0 &&resy(i)==0&&resz(i)==0)
         fprintf(out,'        Fixed end support');
     end
     if(resx(i)==0&&resy(i)==0&&resz(i)==1)
         fprintf(out,'        Pin support');
     end
     if(resx(i)==1&&resy(i)==0&&resz(i)==1)
         fprintf(out,'        Roller support');
     end
     fprintf(out,'\n');
end
fprintf(out,'* NOTE: 1 indicates an unrestrained degree of freedom\n');
fprintf(out,'        0 indicates a  restrained   degree of freedom\n');
fprintf(out,'\nThe vector of restrained degrees of freedom for the structure is\n');
for i = 1 : Model.nBoundary
    fprintf(out,'%d ',Model.bound.ifix(i));
end
fprintf(out,'\n\n');
fprintf(out,'________________________________________________________________________________________________\n');
fprintf(out,'E L E M E N T   I N F O R M A T I O N\n');
fprintf(out,'Element first second Area    Moment of    Young''s    Element            Description\n');
fprintf(out,'        node  node    (A)    Intertia(Iz) modulus(E)  load*    w1       x1       w2       x2\n');
fprintf(out,'------------------------------------------------------------------------------------------------\n');
for e = 1 : Model.nElem
   fprintf(out,'%2d     %2d     %2d     %.3f   %.3e   %.2e     DL  %7.2f %7.2f   %7.2f %7.2f\n',...
       e,Model.geometry.elements(e,1),Model.geometry.elements(e,2),...
       Model.properties.A(e),Model.properties.Iz(e),Model.properties.E(e),...
       wl(e),d1(e),wr(e),d2(e));
end
fprintf(out,'*NOTE: DL = Distributed load w on the element\n');
fprintf(out,'            from a distance of d1 to a distance of d2 from left end\n');
fprintf(out,'________________________________________________________________________________________________\n');
fprintf(out,'DETAILED COMPUTATIONS FOR EACH ELEMENT ARE PRESENTED BELOW\n');
fprintf(out,'{eforce} the vector of equivalent nodal forces for each element along local axes\n');   
fprintf(out,'Element     Vector of Element Forces  {eforce}\n');

for e = 1 : Model.nElem
    fprintf(out,'  %d      ',e);
    for row=1:6
        fprintf(out,'%8.3f   ',eforce(row,e));
    end
    fprintf(out,'\n');
end

end
%--------------------------------------------------------------------------
function printStiff(Solution,fid)
    if nargin < 2
        fid = 1;
    end
    if isfield(Solution,'data')
        printDLine(fid,96);
        printCTitle(fid,'S T I F F N E S S    M A T R I C E S',96);

        nElem = size(Solution.data.ke_local,3);
        for i = 1 : nElem
            printLine(fid,96);
            fprintf(fid,' Element No = %d \n',i);
            printSLine(fid,96);
            fprintf (fid,' - Local Stiffness Matrix: \n');
            disp_mat(Solution.data.ke_local(:,:,i),1,fid,Solution.data.dofList(i,:))
            printSLine(fid,96);

            fprintf (fid,' - Transformation Matrix: \n');
            disp_mat(Solution.data.T6(:,:,i),2,fid,[])        
            printSLine(fid,96);

            fprintf (fid,' - Global Stiffness Matrix: K=T`kT = \n');    
            disp_mat(Solution.data.ke_global(:,:,i),3,fid,Solution.data.dofList(i,:))
        end
        printLine(fid,96);
        printCTitle(fid, 'Assembled Stiffness Matrix', 96)
        printCTitle(fid, '(Before Imposing Boundary Conditions)', 96)
        printSLine(fid,96);
        K0 = Solution.data.K;
        for i = 1 : size(K0,1)
            for j = 1 : size(K0,2)
                fprintf(fid,'%14.3e\t',K0(i,j));
            end
            fprintf(fid,'\n');
        end
        
        printLine(fid,96);
        printCTitle(fid, 'Assembled Stiffness Matrix', 96)
        printCTitle(fid, '(After Imposing Boundary Conditions (Condensed Version))', 96)
        printSLine(fid,96);
        K1 = Solution.data.K_bc;
        for i = 1 : size(K1,1)
            for j = 1 : size(K1,2)
                fprintf(fid,'%14.3e\t',K1(i,j));
            end
            fprintf(fid,'\n');
        end
        
        printDLine(fid,96);
    end 

end
%--------------------------------------------------------------------------
function printResult(Solution,fid)
    
    if nargin < 2
        fid = 1;
    end
    BULLET_CHAR = char(15);
    
    fprintf(fid,'\n\n');
    printDLine(fid,96)
    printCTitle(fid,'O U T P U T     S U M M A R Y',96)
    printDLine(fid,96)
    
    if isfield(Solution,'info')
       fprintf(fid,'%c About: \n',BULLET_CHAR);
       fprintf(fid,'%s\n', Solution.info);
    end
    
    
    uNodal = Solution.nodalDisplacements;
    u = Solution.data.u;
    [nNode,~] = size(uNodal);
    
    printDLine(fid,96)
    printCTitle(fid,'N O D A L   D I S P L A C E M E N T S',96)
    printSLine(fid,96)
    fprintf(fid,'Node        X- displacement           Y- displacement        Z-rotation\n');
    printLine(fid,96)
    for i=1 : nNode
        fprintf(fid,'%2d    %22.10e    %22.10e    %22.10e\n',i,u(3*i-2),u(3*i-1),u(3*i));
    end
    
    lforce = Solution.data.lforce;
    printDLine(fid,96)
    printCTitle(fid,'F O R C E S   A N D   M O M E N T S   I N   E L E M E N T S',96);
    printSLine(fid,96);
    fprintf(fid,'ELEMENT  Axial      Shear     Bending       Axial     Shear     Bending\n');
    fprintf(fid,'         Force      Force     Moment        Force     Force     Moment\n');
    fprintf(fid,'         ----------------------------       ------------------------------\n');
    fprintf(fid,'               at FIRST NODE                       at SECOND NODE\n');
    fprintf(fid,'--------------------------------------------------------------------------\n');
    for e = 1 : size(lforce,2)
        f1= lforce(1,e);
        f2= lforce(2,e);
        f3= lforce(3,e);
        f4= lforce(4,e);
        f5= lforce(5,e);
        f6= lforce(6,e);
        fprintf(fid,'%3d %10.3f %10.3f %10.3f    %10.3f %10.3f %10.3f\n',e,f1,f2,f3,f4,f5,f6);
    end
    
    rf = Solution.reactionForces;
    printDLine(fid,96);
    printCTitle(fid,'S U P P O R T    R E A C T I O N S',96);
    fprintf(fid,'Support node      Rx             Ry          Mz\n');
    printSLine(fid,96);
    for i = 1 : size(rf,1)
        fprintf(fid,'   %2d       %10.3f    %10.3f   %10.3f\n',...
                       rf(i,1),   rf(i,2),  rf(i,3),  rf(i,4));
    end
    
    printDLine(fid,96);
end
%--------------------------------------------------------------------------
function PlanePlot(Model,fig)
    if nargin < 2
        figHandle = figure;
    else
        figHandle = figure(fig);
    end
    clf(figHandle);
    h = subplot(1,1,1);
    
    geometry=Model.geometry;

    coords = geometry.coordinates';
    connect= geometry.elements';
    [ncoord,nnode] = size(coords);
    [nelnode,nelem] = size(connect);
    nelnodes = nelnode*ones(1,nelem);    

    color = 'b';
    fcolor = 'none';
    sw_node = 1;
    sw_elem = 1;
    sw_boun = 1;
    nodeFontSize = 10;
    elemFontSize = 10;
    nodeLabelColor = [1 1 .2];
    elemLabelColor = [1 0 1];
    FACTOR_BC=30;
    
    if isfield (Model,'bound')
        %  Boundary Condition Data
        bc=[Model.bound.resx Model.bound.resy Model.bound.resz]';
        bc = bc(1:2,:);
        nfix = sum(sum(bc==0));
        fixnodes = zeros(2,nfix);
        counter_fix=0;
        for i = 1: nnode
            for j = 1 : 2
                if bc(j,i)==0
                    counter_fix = counter_fix + 1;
                    fixnodes(1,counter_fix)=i;
                    fixnodes(2,counter_fix)=j;
                end
            end
        end
        fixnodes;
    else
       fixnodes = Model.boundary';
    end

    hold(h, 'on')
    axis(h, 'equal')
    grid(h, 'on')

    for lmn = 1:nelem
       for i = 1:nelnodes(lmn)
           x(i,1:2) = coords(1:2,connect(i,lmn));
       end

       pHandle = patch('Parent',h,'Vertices',x,'Faces',[1,2],'FaceColor',fcolor,'EdgeColor',color);

       if sw_elem==1
           lmncoord=coords(:,connect(:,lmn));
           coordp=sum(lmncoord')./nelnodes(lmn);
           text(coordp(1),coordp(2),num2str(lmn),'Parent',h,'BackgroundColor',elemLabelColor,'FontSize',elemFontSize,'FontWeight','bold');
       end    
    end
    if sw_node==1           
       for i=1:nnode
           text(coords(1,i),coords(2,i),num2str(i),'Parent',h,'BackgroundColor',nodeLabelColor,'FontSize',nodeFontSize,'FontWeight','bold');
       end
    end
    %-----------------------------------Boundary
    if sw_boun==1
       xrange=max(coords(1,:))-min(coords(1,:));
       yrange=max(coords(2,:))-min(coords(2,:));
       L1    =max([xrange,yrange]./FACTOR_BC);
       [~,nfix]=size(fixnodes);
       for i=1:nfix
           coordi=coords(:,fixnodes(1,i));
           switch fixnodes(2,i)
               case 1  %x-dir
                   xx=[coordi(1) coordi(1)-L1   coordi(1)-L1  ]; 
                   yy=[coordi(2) coordi(2)+L1/2 coordi(2)-L1/2];
               case 2  %y-dir
                   xx=[coordi(1) coordi(1)+L1/2 coordi(1)-L1/2]; 
                   yy=[coordi(2) coordi(2)-L1   coordi(2)-L1  ];                       
           end
           dt = DelaunayTri(xx',yy');
           triplot(dt,'m','LineWidth',1,'Parent',h);
       end
    end
    %-------------------------------------------
    xlabel(h, '\fontsize{14} X')
    ylabel(h, '\fontsize{14} Y   ', 'rotation', 0);
    %-------------------------------------------
    % expanding limits
    xLim = h.XLim;
    yLim = h.YLim;
    FACTOR_EXP = 0.20;  %20 percent expansion
    
    expansionLengthX = FACTOR_EXP * norm(xLim);
    expansionLengthY = FACTOR_EXP * norm(yLim);
    
    h.XLim = h.XLim + [-1, 1]*expansionLengthX;
    h.YLim = h.YLim + [-1, 1]*expansionLengthY;
    
end
%--------------------------------------------------------------------------
function PlanePlotDeform(Model, Solution, fig)
    if nargin < 3
        figHandle = figure;
    else
        figHandle = figure(fig);
    end
    clf(figHandle);
    h = subplot(1,1,1);
    hold(h, 'on')
    axis(h, 'equal')
    grid(h, 'on')
    
    geometry=Model.geometry;

    coords = geometry.coordinates';
    connect= geometry.elements';
    [ncoord,nnode] = size(coords);
    [nelnode,nelem] = size(connect);
    nelnodes = nelnode*ones(1,nelem);    

    color = 'b';
    fcolor = 'none';
    sw_node = 1;
    sw_elem = 1;
    sw_boun = 1;
    nodeFontSize = 10;
    elemFontSize = 10;
    nodeLabelColor = [1 1 .2];
    elemLabelColor = [1 0 1];
    FACTOR_BC=30;
    
    if isfield (Model,'bound')
        %  Boundary Condition Data
        bc=[Model.bound.resx Model.bound.resy Model.bound.resz]';
        bc = bc(1:2,:);
        nfix = sum(sum(bc==0));
        fixnodes = zeros(2,nfix);
        counter_fix=0;
        for i = 1: nnode
            for j = 1 : 2
                if bc(j,i)==0
                    counter_fix = counter_fix + 1;
                    fixnodes(1,counter_fix)=i;
                    fixnodes(2,counter_fix)=j;
                end
            end
        end
        fixnodes;
    else
       fixnodes = Model.boundary';
    end

    hold(h, 'on')
    axis(h, 'equal')
    grid(h, 'on')

    for lmn = 1:nelem
       for i = 1:nelnodes(lmn)
           x(i,1:2) = coords(1:2,connect(i,lmn));
       end

       pHandle = patch('Parent',h,'Vertices',x,'Faces',[1,2],'FaceColor',fcolor,'EdgeColor',color,...
                       'LineWidth',1.5);

       if sw_elem==1
           lmncoord=coords(:,connect(:,lmn));
           coordp=sum(lmncoord')./nelnodes(lmn);
           text(coordp(1),coordp(2),num2str(lmn),'Parent',h,'BackgroundColor',elemLabelColor,'FontSize',elemFontSize,'FontWeight','bold');
       end    
    end
    if sw_node==1           
       for i=1:nnode
           text(coords(1,i),coords(2,i),num2str(i),'Parent',h,'BackgroundColor',nodeLabelColor,'FontSize',nodeFontSize,'FontWeight','bold');
       end
    end
    %-----------------------------------Boundary
    if sw_boun==1
       xrange=max(coords(1,:))-min(coords(1,:));
       yrange=max(coords(2,:))-min(coords(2,:));
       L1    =max([xrange,yrange]./FACTOR_BC);
       [~,nfix]=size(fixnodes);
       for i=1:nfix
           coordi=coords(:,fixnodes(1,i));
           switch fixnodes(2,i)
               case 1  %x-dir
                   xx=[coordi(1) coordi(1)-L1   coordi(1)-L1  ]; 
                   yy=[coordi(2) coordi(2)+L1/2 coordi(2)-L1/2];
               case 2  %y-dir
                   xx=[coordi(1) coordi(1)+L1/2 coordi(1)-L1/2]; 
                   yy=[coordi(2) coordi(2)-L1   coordi(2)-L1  ];                       
           end
           dt = DelaunayTri(xx',yy');
           triplot(dt,'m','LineWidth',1,'Parent',h);
       end
    end
    %-------------------------------------------
    xlabel(h, '\fontsize{14} X')
    ylabel(h, '\fontsize{14} Y   ', 'rotation', 0);
    %-------------------------------------------
    % Plotting deformations
    defColor = 'r';
    u = Solution.nodalDisplacements;
    u = u(:,1:Model.nDim);
    maxU = max(max(abs(u)));
    
    % exaggeration factors:
    % expanding limits
    xLim = h.XLim;
    yLim = h.YLim;
    FACTOR_MAGNIFYING = 0.10;  %10 percent of axis limits
    
    exaggFactor = FACTOR_MAGNIFYING * max(norm(xLim),norm(yLim)) / maxU;
    
    coords = (geometry.coordinates + exaggFactor.*u)';
    
    for lmn = 1:nelem
       for i = 1:nelnodes(lmn)
           x(i,1:2) = coords(1:2,connect(i,lmn));
       end

       pHandle = patch('Parent',h,'Vertices',x,'Faces',[1,2],'FaceColor',fcolor,'EdgeColor',defColor,...
                       'LineStyle','--');  
    end
     title(h,{'\fontsize{15}\bf Deformed Shape '; ...
              ['\color{red}\fontsize{11}\rm\it Automatic Magnifying Factor = ',num2str(round(exaggFactor))]});
%            ['\color[rgb]{0 .5 .5} \fontsize{11}\rm\it Magnifying Factor = ',num2str(exaggFactor)]});
    %-------------------------------------------
    % adjusting (expanding) limits
    xLim = h.XLim;
    yLim = h.YLim;
    FACTOR_EXP = 0.20;  %20 percent expansion
    
    expansionLengthX = FACTOR_EXP * norm(xLim);
    expansionLengthY = FACTOR_EXP * norm(yLim);
    
    h.XLim = h.XLim + [-1, 1]*expansionLengthX;
    h.YLim = h.YLim + [-1, 1]*expansionLengthY;
end
%--------------------------------------------------------------------------
function printLine(fid,n)
    if nargin < 2
        n = 80; 
    end
    if nargin < 1
        fid = 1;
    end
    fprintf(fid,'%s\n',repmat('_',1,n));
end
function printDLine(fid,n)
    if nargin < 2
        n = 80; 
    end
    if nargin < 1
        fid = 1;
    end
    fprintf(fid,'%s\n',repmat('=',1,n));
end
function printSLine(fid,n)
    if nargin < 2
        n = 80; 
    end
    if nargin < 1
        fid = 1;
    end
    fprintf(fid,'%s\n',repmat('-',1,n));
end
function printCTitle(fid,text,n)
    if nargin < 3
        n = 80;
    end
    if nargin < 2
        text =''; 
    end
    if nargin < 1
        fid = 1;
    end
    tLen = length(text);
    space = floor((n-tLen)/2);
    fprintf(fid,'%s\n',[repmat(' ',1,space),text]);
end
function disp_mat(mat,mode,fid,dofList)
format compact

if mode==1 
    format shortEng
    
    fprintf(fid,'     +----                                                    ----+  \n');
    fprintf(fid,'     | EA/L    0          0         -EA/L      0          0       |  \n');
    fprintf(fid,'     | 0       12EI/L^3   6EI/L^2    0        -12EI/L^3   6EI/L^2 |  \n');
    fprintf(fid,'     | 0       6EI/L^2    4EI/L      0        -6EI /L^2   2EI/L   |  \n');
    fprintf(fid,'     |-EA/L    0          0          EA/L      0          0       | =\n');
    fprintf(fid,'     | 0      -12EI/L^3  -6EI/L^2    0         12EI/L^3  -6EI/L2  |  \n');
    fprintf(fid,'     | 0       6EI/L^2    2EI/L      0        -6EI /L^2   4EI/L   |  \n');
    fprintf(fid,'     +----                                                    ----+  \n');
    for i = 1 : 6
        for j = 1 : 6
            fprintf(fid,'%14.3e\t',mat(i,j));
        end
        fprintf(fid,'\n');
    end

elseif mode==2
    format short g
    c=mat(1,1);
    s=mat(1,2);
    fprintf(fid,'     +---                  ---+   +---                               ---+\n');
    fprintf(fid,'     |  c   s   0   0   0   0 |   | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',c ,s ,0,0,0,0);
    fprintf(fid,'     | -s   c   0   0   0   0 |   | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',-s,c ,0,0,0,0);
    fprintf(fid,'     |  0   0   1   0   0   0 |   | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',0 ,0 ,1,0,0,0);
    fprintf(fid,'     |  0   0   0   c   s   0 | = | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',0 ,0 ,0,c,s,0);
    fprintf(fid,'     |  0   0   0  -s   c   0 |   | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',0 ,0 ,0,-s,c,0);
    fprintf(fid,'     |  0   0   0   0   0   1 |   | %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f |\n',0 ,0 ,0,0,0,1);
    fprintf(fid,'     +---                  ---+   +---                               ---+\n');

elseif mode==3
    fprintf(fid,'     +----                                                                  ----+ \n');
    for i = 1 : 6
        fprintf(fid, '     |');
        for j = 1 : 6
            fprintf(fid, ' k(%2d,%2d) \t',dofList(i),dofList(j));
        end
        % fprintf(fid, '| \n');
        if i == floor(size(mat,1)/2)+1
            fprintf(fid, '|  =  \n');
        else
            fprintf(fid, '| \n');
        end
    end
    fprintf(fid,'     +----                                                                  ----+ \n');
    
    for i = 1 : 6
        for j = 1 : 6
            fprintf(fid,'%14.3e\t',mat(i,j));
        end
        fprintf(fid,'\n');
    end
end


fprintf('\n')
format long 
end
%--------------------------------------------------------------------------