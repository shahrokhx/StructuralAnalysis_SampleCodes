%%                P L A N A R     T R U S S     S O L V E R                 
%__________________________________________________________________________

%               Developed by  SHAHROKH SHAHI (www.sshahi.com) 
%                   Georgia Institute of Technology
%__________________________________________________________________________
%
%
%  Usage: 
%  The only inputs that you need to provide is:
%  (1) input file name (LINE 28) 
%  (2) output file name (LINE 29)
%  (3) a (0/1) switch value for displaying stiffness matrices (LINE 32) 
%  

%% Inputs
function Solution = planarTrussSolver(inpFileName)

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
figureNumber = 1;
 PlanePlot(Model,figureNumber);

% analysing the model
Solution = trussAnalysis(Model);

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

end

%% Functions
%--------------------------------------------------------------------------
function Solution = trussAnalysis(Model)
% PURPOSE: Analysing trusses with explicit stiffness matrices    
%
% INPUT(S): 
%   - Model: a structure including the model data
%
% OUTPUT(S):
%   - Solution: a structure including the solution of the analysis

% storing information in Solution structs
time = clock;
Solution.info = sprintf(['Solution obtained on %s %d:%d:%d',...
                         '\nExplicit Stiffness Matrices' ],...
                         date,time(4),time(5),floor(time(6)));
% Pre-processing
%--------------------------------------------------------------------------
%                        I N P U T     D A T A
%--------------------------------------------------------------------------

% control variables

nNode   = Model.nNode;
nElem   = Model.nElem;
nDof    = Model.nDof;

%nodal coordinates
coordinates = Model.geometry.coordinates;

% nodal forces
F = zeros(size(coordinates));
for i = 1 : Model.nNode
    node = Model.loading.nodalForces(i,1);
    F(node, 1:Model.nDim) = Model.loading.nodalForces(i,2:end);
end

% reshape to vector form
Feq=reshape(F',1,nDof*nNode)';  

% elements connectivity
elements = Model.geometry.elements;

% section/materaial properties
A = zeros(1, Model.nElem);
E = zeros(1, Model.nElem);
for i = 1 : Model.nElem
    E(i) = Model.sections(Model.elemSectId(i), 1);
    A(i) = Model.sections(Model.elemSectId(i), 2);
end

% boundary condition
boundary = Model.boundary';

%--------------------------------------------------------------------------
%                         P R O C E S S I N G
%--------------------------------------------------------------------------
% initializing the loop variables
K = zeros(nDof*nNode);
k(:,:,nElem) = zeros(2 * nDof);
dofList      = zeros(nElem,4); % to store the assembly indices
L = zeros(1,nElem);
k11 = zeros(1,nElem);
% loop over elements
for i = 1 : nElem
    a = coordinates(elements(i,1),1:2); % start node coordinates
    b = coordinates(elements(i,2),1:2); % end node coordinates
    
    % calculating length
    L(i) = sqrt((a(1)-b(1)).^2 + (a(2)-b(2)).^2);  
    
    % sin and cos
    c = (b(1) - a(1))/L(i);  % cos(t) = (x2 - x1)/L
    s = (b(2) - a(2))/L(i);  % sin(t) = (y2 - y1)/L
    
    % TODO: this only works for 2D
    lambda = [c*c c*s; c*s s*s];
    T = [lambda -lambda; -lambda lambda];
    
    k11(i)   = E(i)*A(i)/L(i);
    k(:,:,i) = E(i)*A(i)/L(i)   *    T;
    
    % TODO : this only works for 2D
    T_sigma(:,:,i) = [c s 0 0; 0 0 c s];
    
    % node i --------> dof ndof*i - ndof + 1 : ndof*i
    % e.g.  node 1 --> dof 1 : 2
    %       node 3 --> dof 5 : 6
    inode = elements(i,1);
    jnode = elements(i,2);
    
    index = [nDof*(inode-1)+1:nDof*inode,nDof*(jnode-1)+1:nDof*jnode];
    dofList(i,:) = index;
    
    K(index,index) = K(index,index) + k(:,:,i);
    
end
Solution.data.elementStiff = k;
Solution.data.globalStiff = K;
Solution.data.globalForce = Feq;
Solution.data.dofList = dofList;
Solution.data.rotation = T_sigma;
Solution.data.k11 = k11;

% imposing displacement boundary conditions
debc = double(nDof) .* (boundary(1,:) - 1) + boundary(2,:);
ebcVals = boundary(3,:)';
[dofs, rf] = solution(K, Feq, debc, ebcVals);
Solution.nodalDisplacements = reshape(dofs,Model.nDim,Model.nNode)';

Solution.reactionForces = zeros(Model.nNode, Model.nDim);
for i = 1 : size(boundary,2)
    node = boundary(1,i);
    dof = boundary(2,i);
    Solution.reactionForces(node,dof) = rf(i);
end


df = setdiff(1:length(Feq), debc);
Solution.data.K_bc = K(df,df);
%--------------------------------------------------------------------------
%                      P O S T - P R O C E S S I N G
%--------------------------------------------------------------------------    
% calculating element stresses

Solution.elementStresses.sigma11 = zeros(nElem,1);
Solution.internalForces.axial = zeros(nElem,1);

for i = 1 : nElem    
    
    inode = elements(i,1);
    jnode = elements(i,2);
    
    index = [nDof*(inode-1)+1:nDof*inode,nDof*(jnode-1)+1:nDof*jnode];
    
    sigma = E(i) * [-1/L(i), 1/L(i)] * T_sigma(:,:,i) * dofs(index);
    
    Solution.elementStresses.sigma11(i) = sigma;
    Solution.internalForces.axial(i,1) = sigma * A(i);
end

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
function [d, rf] = solution(K, R, debc, ebcVals)
    dof = length(R);
    df = setdiff(1:dof, debc);
    Kf = K(df, df);
    Rf = R(df) - K(df, debc)*ebcVals;
    dfVals = Kf\Rf;
    d = zeros(dof,1);
    d(debc) = ebcVals;
    d(df) = dfVals;
    rf = K(debc,:)*d - R(debc);
end
%--------------------------------------------------------------------------
function Model = inpFileReader(inpFileName)
inp   = fopen(inpFileName,'r');
%--------------------------------------------------------------------------
%       R E A D I N G     G E N E R A L      I N F O R M A T I O N
%--------------------------------------------------------------------------
mat             = fscanf(inp,'%f',2); 
Model.nNode     = mat(1);              
Model.nElem     = mat(2);             
Model.nDim      = 2;
Model.nDof      = 2;
Model.nElemNode = 2;
mat             = (fscanf(inp,'%d %f %f %d %d  %f %f',[7,Model.nNode]))';
Xval            = mat(:,2);
Yval            = mat(:,3);
resx            = mat(:,4);
resy            = mat(:,5);
fx              = mat(:,6);
fy              = mat(:,7);
mat             = (fscanf(inp,'%d %d %d %f %f',[5,Model.nElem]))';
FirstNode       = mat(:,2);
SecondNode      = mat(:,3);
E               = mat(:,4);
A               = mat(:,5);

F = zeros(Model.nNode*Model.nDof,1); 
for i = 1 : Model.nNode
  F(2*i-1) = fx(i);  
  F(2*i  ) = fy(i); 
end

% 0 = RESTRAINED degree of freedom
% 1 = UNRESTRAINED degree of freedom
bound = [];
for i = 1 : Model.nNode
    if resx(i)==0
        bound = [bound; [i, 1, 0.0]];
    end 
    if resy(i)==0
        bound = [bound; [i, 2, 0.0]];
    end    
end
nBoundary = size(bound,1);


Model.geometry.coordinates   = [Xval, Yval];
Model.geometry.elements      = [FirstNode, SecondNode];
Model.nNodalForce            = Model.nNode;
Model.nTractionForce         = 0;
Model.nBodyForce             = 0;
Model.loading.nodalForces    = [(1:Model.nNode)', fx, fy];
Model.loading.tractionForces = [];
Model.loading.bodyForces     = [];
Model.nSection               = Model.nElem;
Model.nMaterialProp          = 2;
Model.elemSectId             = 1 : Model.nElem;
Model.sections               = [E, A];
Model.nBoundary              = nBoundary;
Model.boundary               = bound;

end
%--------------------------------------------------------------------------
function printSummary(Model,fid)
% PURPOSE: reads input files and generates a structure including all data    

    if nargin < 2 
        fid = 1;
    end
    BULLET_CHAR = char(15);
    
    printDLine(fid)
    tText = 'A    S U M M A R Y    O F    T H E    I N P U T    M O D E L';
    printCTitle(fid,tText)
    printDLine(fid)


    if isfield(Model,'info')
       fprintf(fid,'%c About: \n',BULLET_CHAR);
       fprintf(fid,'%s\n', Model.info);
    end
    if isfield(Model,'analysisType')
       fprintf(fid,'\n%c Structure Type: \n',BULLET_CHAR);
       fprintf(fid,'%s\n', Model.analysisType);
    end
    printLine(fid)

    
    fprintf(fid,'%c Control Variables: \n',BULLET_CHAR);
    fprintf(fid,'   - Dimension:             %dD\n' ,Model.nDim);
    fprintf(fid,'   - Number of DOF(s):      %d \n' ,Model.nDof);
    fprintf(fid,'\n');
    fprintf(fid,'   - Number of Nodes:    %4d \n' ,Model.nNode);
    fprintf(fid,'   - Number of Elements: %4d \n' ,Model.nElem);
    fprintf(fid,'   - Nodes per Elements:    %d\n',Model.nElemNode);
    fprintf(fid,'\n');
    fprintf(fid,'   - Number of Sections:    %d \n',Model.nSection);
    printDLine(fid)

    
    printCTitle(fid,'N O D A L    I N F O R M A T I O N')
    printLine(fid)
    headerTxt = ['Node     Coordinates          Nodal Loads       ',...
                 '         Nodal Restraints    \n'];
    fprintf(fid,headerTxt);
    printSLine(fid)
    switch Model.nDim
        case 1 
            headerTxt1 =' ID          X               ';
            headerTxt2 =' Fx                      ';
            headerTxt3 =' Rx          \n';
        
        case 2
            headerTxt1 =' ID       X     Y            ';
            headerTxt2 =' Fx       Fy            ';
            headerTxt3 =' Rx   Ry     \n';
            
        case 3
            headerTxt1 =' ID     X    Y    Z          ';
            headerTxt2 =' Fx       Fy       Fz   ';
            headerTxt3 =' Rx   Ry   Rz \n';
    end
    fprintf(fid,[headerTxt1,headerTxt2,headerTxt3]);
    printLine(fid)
    
    % nodal forces to displaying
    F = zeros(Model.nNode,Model.nDof);
    for i = 1 : Model.nNodalForce
        node = Model.loading.nodalForces(i,1);
        F(node, 1:Model.nDof) = Model.loading.nodalForces(i,2:end);
    end
        
    % printing table contents
    for i = 1 : Model.nNode
        coord = Model.geometry.coordinates(i,:);
        fprintf(fid,'%2d  ',i);
        for j = 1 : length(coord)
            fprintf(fid,'%7.2f  ',coord(j));
        end
        fprintf(fid,'\t');
        for j = 1 : size(F,2)
            fprintf(fid,'%10.2f  ',F(i,j));
        end
        
        fprintf(fid,'\t');
        index = find(Model.boundary(:,1)==i);
        dx='-';dy='-';dz='-';bc='';
        dx='1';dy='1';dz='1';bc='';
        
        if ~isempty(index)
            for k = 1 : length(index)
                if Model.boundary(index(k),2)==1
                    dx=num2str(Model.boundary(index(k),3));
                elseif Model.boundary(index(k),2)==2
                    dy=num2str(Model.boundary(index(k),3));
                elseif Model.boundary(index(k),2)==3
                    dz=num2str(Model.boundary(index(k),3));
                end
             end
             if length(index)==1
                 bc = 'Roller ';
             elseif length(index)==2
                 bc = 'Pin';
             elseif length(index)==3
                 bc = 'Fixed';
             end
        end
        presc = [dx, dy, dz]; %prescribed values
        for j = 1 : Model.nDof
            fprintf(fid,'   %s',presc(j));
        end
        fprintf(fid,'       %s\n',bc);
    end
    printDLine(fid);
    
    printCTitle(fid,'E L E M E N T S     I N F O R M A T I O N')
    printLine(fid)
    fprintf(fid,' ID\t ');
    for i = 1 : Model.nElemNode
        fprintf(fid,'Node%d\t',i);
    end
    fprintf(fid,'     Section No.     Material Properties\n');
    printLine(fid);
    
    for i = 1 : Model.nElem
        fprintf(fid,'%2d\t ',i);
        for j = 1 : Model.nElemNode
            fprintf(fid,'%2d\t\t',Model.geometry.elements(i,j));
        end
        fprintf(fid,'\t');
        fprintf(fid,'\t%2d\t\t\t[',Model.elemSectId(i));
        for j = 1 : Model.nMaterialProp
            fprintf(fid,'%g\t\t',Model.sections(Model.elemSectId(i),j));
        end
        fprintf(fid,']\n');
    end
    
    printDLine(fid);
    
end
%--------------------------------------------------------------------------
function printStiff(Solution,fid)

    if nargin < 2
        fid = 1;
    end
    
    if isfield(Solution,'data')
        fprintf(fid,'\n\n\n\n');
        % element stiffness
        printDLine(fid);
        if isfield(Solution.data,'elementStiff')
            printCTitle(fid,'S T I F F N E S S     M A T R I C E S')
            
    
            [n,m,nElem] = size(Solution.data.elementStiff);
            for iElem = 1 : nElem
                stif = Solution.data.elementStiff(:,:,iElem);
                index= Solution.data.dofList(iElem,:);
                rotMat = Solution.data.rotation(:,:,iElem);
                c = rotMat(1,1);
                s = rotMat(1,2);
                k = Solution.data.k11(iElem);
                
                printLine(fid);
                fprintf(fid, 'Element No. [%2d]\n',iElem);
                printSLine(fid);
                fprintf (fid,' - Local Stiffness Matrix: \n');
                fprintf(fid,'     +----                ----+   \n');
                fprintf(fid,'     | EA/L   0    -EA/L    0 |   \n');
                fprintf(fid,'     | 0      0     0       0 |   \n');
                fprintf(fid,'     |-EA/L   0     EA/L    0 | = \n');
                fprintf(fid,'     | 0      0     0       0 |   \n');
                fprintf(fid,'     +----                ----+   \n');
                stifLocal = zeros(4);
                stifLocal([1,3],[1,3]) = [k,-k;-k,k];
                for i = 1 : n
                    for j = 1 : m
                        fprintf(fid,'%14.3e\t',stifLocal(i,j));
                    end
                    fprintf(fid,'\n');
                end
                
                printSLine(fid);
                fprintf (fid,' - Transformation Matrix: \n');
                fprintf(fid,'     +---           ---+   +---                    --+\n');
                fprintf(fid,'     |  c   s   0   0  |   | %5.2f %5.2f %5.2f %5.2f |\n',c ,s , 0,0);
                fprintf(fid,'     | -s   c   0   0  |   | %5.2f %5.2f %5.2f %5.2f |\n',-s,c , 0,0);
                fprintf(fid,'     |  0   0   c   s  | = | %5.2f %5.2f %5.2f %5.2f |\n',0 ,0 , c,s);
                fprintf(fid,'     |  0   0  -s   c  |   | %5.2f %5.2f %5.2f %5.2f |\n',0 ,0 ,-s,c);
                fprintf(fid,'     +---           ---+   +---                   ---+\n');
                
                printSLine(fid);
                fprintf (fid,' - Global Stiffness Matrix: K=T`kT = \n');
                fprintf(fid,'     +----                                    ----+ \n');
    
                for i = 1 :size(stif,1)
                    fprintf(fid, '     |');
                    for j = 1 : size(stif,2)
                        fprintf(fid, ' k(%2d,%2d)  ',index(i),index(j));
                    end
                    if i == floor(size(stif,1)/2)+1
                        fprintf(fid, '|  =  \n');
                    else
                        fprintf(fid, '| \n');
                    end
                end
                fprintf(fid,'     +----                                    ----+ \n\n');


                for i = 1 : n
                    for j = 1 : m
                        fprintf(fid,'%14.3e\t',stif(i,j));
                    end
                    fprintf(fid,'\n');
                end
            end
        end
        printLine(fid);
        
        % global assembled stiffness
        if isfield(Solution.data,'globalStiff')
            printCTitle(fid, 'Assembled Stiffness Matrix');
            printCTitle(fid, '(Before Imposing Boundary Conditions)');
            printSLine(fid);
            stif = Solution.data.globalStiff;
            [n,m] = size(stif);
            for i = 1 : n
                for j = 1 : m
                    fprintf(fid,'%14.3e\t',stif(i,j));
                end
                fprintf(fid,'\n');
            end
        end
        printLine(fid);
        
        if isfield(Solution.data,'K_bc')
            printCTitle(fid, 'Assembled Stiffness Matrix');
            printCTitle(fid, '(After Imposing Boundary Conditions (Condensed Version))');
            printSLine(fid);
            stif = Solution.data.K_bc;
            [n,m] = size(stif);
            for i = 1 : n
                for j = 1 : m
                    fprintf(fid,'%14.3e\t',stif(i,j));
                end
                fprintf(fid,'\n');
            end
        end
%         printLine(fid);
%         
%         
%         % global forces
%         if isfield(Solution.data,'globalForce')
%             fprintf(fid, 'Global Nodal Forces\n');
%             F = Solution.data.globalForce;
%             [n,m] = size(F);
%             for i = 1 : n
%                 fprintf(fid,'[%3d] \t',i);
%                 for j = 1 : m
%                     fprintf(fid,'%14.3e\t',F(i,j));
%                 end
%                 fprintf(fid,'\n');
%             end
%         end
        
        printDLine(fid);
    end
end
%--------------------------------------------------------------------------
function printResult(Solution,fid)
    
    if nargin < 2
        fid = 1;
    end
    
    BULLET_CHAR = char(15);
    
    fprintf(fid,'\n\n\n\n');
    printDLine(fid)
    printCTitle(fid,'O U T P U T     S U M M A R Y')
    printDLine(fid)
    
    if isfield(Solution,'info')
       fprintf(fid,'%c About: \n',BULLET_CHAR);
       fprintf(fid,'%s\n', Solution.info);
    end
    printLine(fid)
    
    
    if isfield(Solution, 'nodalDisplacements')
        u = Solution.nodalDisplacements;
        [nNode,nDof] = size(u);
        printCTitle(fid,'Nodal Displacements');
        printSLine(fid);
        for i = 1 : nNode
            formatStr = repmat('%16.8f\t',1,nDof);
            fprintf(fid,['[%3d]   ',formatStr,'\n'],i,u(i,:));
        end
    end
    printLine(fid)
    if isfield(Solution, 'reactionForces')
        rf = Solution.reactionForces;
        [nNode,nDof] = size(rf);
        printCTitle(fid,'Reaction Forces');
        printSLine(fid);
        for i = 1 : nNode
            formatStr = repmat('%10.5f\t',1,nDof);
            fprintf(fid,['[%3d]   ',formatStr,'\n'],i,rf(i,:));
        end
    end
    printLine(fid)
    if isfield(Solution, 'elementStresses')
        stress = Solution.elementStresses.sigma11;
        [nElem,nComponenet] = size(stress);
        printCTitle(fid,'Stresses');
        printSLine(fid);
        for i = 1 : nElem
            formatStr = repmat('%10.5f\t',1,nComponenet);
            fprintf(fid,['[%3d]   ',formatStr,'\n'],i,stress(i,:));
        end
    end
    printLine(fid)
    if isfield(Solution, 'internalForces')
        internalForces = Solution.internalForces.axial;
        [nElem,nComponenet] = size(internalForces);
        printCTitle(fid,'Internal Forces (Axial)');
        printSLine(fid);
        for i = 1 : nElem
            formatStr = repmat('%10.5f\t',1,nComponenet);
            fprintf(fid,['[%3d]   ',formatStr,'\n'],i,internalForces(i,:));
        end
    end
    printDLine(fid);
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
%--------------------------------------------------------------------------