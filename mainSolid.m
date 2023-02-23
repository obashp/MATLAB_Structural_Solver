clc
clear all

wrkDir = './' ;
problemString = 'rigidV3' ;
elemType = '4Quad' ;
problemType = '2D' ;

% Coordinate, connectivity and boundary data
% dataStr = strcat(wrkDir,'Data_flexible_filament.mat') ;
% load(dataStr);
% ndof = size(crd,1) ;
% conn = BCStructure;
% crdSolid = crd(unique(reshape(BCStructure,[],1)),2:end);


nElem_x = 64 ;      % Number of elements in X
nElem_y = 16 ;      % Number of elements in Y
nen   = 4 ;         % Number of nodes per element
nElem = nElem_x * nElem_y ;
xmax = 0.35;
ymax = 0.02;
dx = xmax/nElem_x; dy = ymax/nElem_y;

x0 = 0.2 + 0.05*cos(asin(0.2));
y0 = 0.2 ;

[xcrd, ycrd] = ndgrid(x0:dx:x0+xmax,y0-0.01:dy:y0+0.01); 
crd = [xcrd(:), ycrd(:)];
ndof = size(crd,1);

targcrd = [x0+xmax,y0];
xtarg = find(abs(xcrd - targcrd(1)) < 1e-12);
ytarg = find(abs(ycrd - targcrd(2)) < 1e-12);
targ = intersect(xtarg,ytarg);

cnn = zeros((nElem_x)*(nElem_y),nen);
for i=1:nElem_y
    for j=1:nElem_x
        ielem = (i-1)*(nElem_x)+j;
        inode = (i-1)*(nElem_x+1)+j;
        % Connectivity
        cnn(ielem,:) = [inode inode+1 inode+(nElem_x+2) inode+(nElem_x+1)];
    end
end

Bc_y0 = [2:nElem_x+1]';
Bc_xL = [2*(nElem_x+1):nElem_x+1:(nElem_y+1)*(nElem_x+1)]';
Bc_yL = [nElem_y*(nElem_x+1)+nElem_x:-1:nElem_y*(nElem_x+1)+2]';
Bc_x0 = [(nElem_y)*(nElem_x+1)+1:-(nElem_x+1):nElem_x+2]';
Bc_x0(end+1) = 1;

theta = atan((crd(Bc_x0,2) - 0.2)./(crd(Bc_x0,1)-0.2));
crd(Bc_x0,1) = 0.2 + 0.05*cos(theta);
deltax = max(crd(Bc_x0,1))- min(crd(Bc_x0,1));
crd(setdiff(1:length(crd(:,1)),Bc_x0),1) = crd(setdiff(1:length(crd(:,1)),Bc_x0),1)+deltax;

StoreVTK('Solid_Stress','2D',ndof,nElem,nen,crd,cnn);

% Nonlinear iteration data:
solver.nLIterMin = 5 ;
solver.nLIterMax = 10 ;
solver.dt = 1e-5 ;
solver.maxSteps = 75000;
solver.rhoinfty = 1.0 ;
solver.nLTol = 1e-4 ;
solver.outFreq = 100 ;
solver.intgOutFreq = 1 ;

% Solid properties
solid.dens = 1000.0 ;
solid.nu   = 0.4 ;
solid.gravFrc = [0,-2.0];
solid.E    = 1.4e6;
solid.mu   = 0.5*solid.E/(1+solid.nu);
solid.lambda = solid.nu*solid.E/((1+solid.nu)*(1-2*solid.nu));

% Set boundary conditions, all surfaces are stress free, and beam is
% cantilevered at one end.
solid.Bc_LeftN = size(unique(Bc_x0),1);
solid.Bc_LeftNs = unique(Bc_x0);
solid.Bc_LeftV = [0 0];

solid.DirichletU = [solid.Bc_LeftNs];
solid.DirichletUval = [solid.Bc_LeftV(1).*ones(solid.Bc_LeftN,1)];
solid.DirichletV = [solid.Bc_LeftNs];
solid.DirichletVval = [solid.Bc_LeftV(2).*ones(solid.Bc_LeftN,1)];

% Gen-alpha time integration parameters
pmc.alphaM = 0.5*(3.0-solver.rhoinfty)/(1+solver.rhoinfty) ;
pmc.alpha = 1.0/(1+solver.rhoinfty) ;
pmc.gamma2 = 0.5 + pmc.alphaM - pmc.alpha ;
pmc.gamma1 = 0.5 + 0.5*(pmc.alphaM - pmc.alpha);
% pmc.beta  = 0.25*(1+pmc.alphaM-pmc.alpha)^2;
% pmc.alphaM = 1.0;
% pmc.alpha = 1.0;
% pmc.gamma1 = 1.0;
% pmc.gamma2 = 1.0;
% % pmc.beta = 0.5;

% Solid variables
Sol.du = zeros(ndof,2,1);
Sol.u = zeros(ndof,2,1);        %Displacement vector (u_n+1)
Sol.v = zeros(ndof,2,1);        %Velocity vector (v_n+1)
Sol.uDot = zeros(ndof,2,1);     %time derivative of displacement vector
Sol.vDot = zeros(ndof,2,1);     %Acceleration of point (vDot_n+1)



Sol.uAlpha = zeros(ndof,2,1) ;  % u-alpha
Sol.vAlpha = zeros(ndof,2,1) ;  % v-alpha
Sol.vDotAlpha = zeros(ndof,2,1) ;%vDot-alpha
Sol.uPrev = Sol.u ;             % u_n
Sol.vPrev = Sol.v ;             % v_n
Sol.vDotPrev = Sol.vDot ;       % vDot_n

filename1 = sprintf('%s/%s.oisd',wrkDir,problemString);
fileId1 = fopen(filename1,'w');

Sol.crd = crd;

time = 1:solver.maxSteps;
disp = zeros(solver.maxSteps,2);
timeStep = 1;

for timeStep = 1:solver.maxSteps
    fprintf('Time step:%d\n',timeStep);
    
%     Predict the solution
%     Sol.vDot = 0 ;  
%     Sol.v = Sol.vPrev + solver.dt*(1.0-pmc.gamma)*pmc.alpha*Sol.vDotPrev;
%     Sol.u = Sol.uPrev + solver.dt*(1.0-pmc.gamma)*pmc.alpha*Sol.vPrev;
    Sol.vDot = Sol.vDotPrev;
    Sol.v = Sol.vPrev + Sol.vDotPrev*solver.dt;
    Sol.u = Sol.uPrev + Sol.vPrev*solver.dt;

%     Sol.u = Sol.uPrev + Sol.vPrev*solver.dt + solver.dt^2 * (0.5 - pmc.beta/pmc.gamma) * Sol.vDotPrev;
%     Sol.u = Sol.uPrev;
%     Sol.v = Sol.vPrev;
%     Sol.vDot = (1.0 - 1.0/pmc.gamma)*Sol.vDotPrev;

    % Interpolate for alpha values for Gen-alpha
%     Sol.uAlpha = Sol.u + pmc.alpha.*(-Sol.u + Sol.uPrev) ;
%     Sol.vAlpha = Sol.v + pmc.alpha.*(-Sol.v + Sol.vPrev) ;
%     Sol.vDotAlpha = Sol.vDot + pmc.alphaM.*(-Sol.vDot + Sol.vDotPrev) ;

    Sol.uAlpha = (1.0-pmc.alpha)*Sol.uPrev + pmc.alpha*Sol.u;
    Sol.vDotAlpha = (1.0-pmc.alphaM)*Sol.vDotPrev + pmc.alpha*Sol.vDot;
                                          
    % Nonlinear iterations start here
    for nLiter = 1:solver.nLIterMax    
        % Solve Elasticity equation        
        [Sol, HEnormIndicator] = linearElastic(solver, solid, pmc, Sol, cnn, crd, ... 
                                              elemType, ndof, nen, nElem);
        
        
        % Check convergence criteria
        if (HEnormIndicator < solver.nLTol)
            break;
        end
    end
    
    % Update the solution
%     Sol.u = Sol.uPrev + (1/(1.0-pmc.alpha))*( Sol.uAlpha - Sol.uPrev );
%     Sol.vDot = Sol.vDotPrev + (1/(1.0-pmc.alphaM))*( Sol.vDotAlpha - Sol.vDotPrev ) ;

%     Sol.v = Sol.vPrev + (pmc.gamma*Sol.vDot + (1.0-pmc.gamma)*Sol.vDotPrev);
    Sol.u = (1.0/pmc.alpha)*Sol.uAlpha + (1.0 - 1.0/pmc.alpha)*Sol.uPrev;
    Sol.vDot = (1.0/pmc.alphaM)*Sol.vDotAlpha + (1.0 - 1.0/pmc.alphaM)*Sol.vPrev;
    Sol.v = (1.0/pmc.alpha)*Sol.vAlpha + (1.0 - 1.0/pmc.alpha)*Sol.vPrev;

%     
%     
%     %Print displacement of target node
    time(timeStep) = timeStep*solver.dt;
    disp(timeStep,1) = Sol.u(targ,1);
    disp(timeStep,2) = Sol.u(targ,2);
    vel(timeStep,1) = Sol.v(targ,1);
    vel(timeStep,2) = Sol.v(targ,2);
    fprintf('displacements = %e %e, \n', Sol.u(targ,1),Sol.u(targ,2));
     
    % Copy current variables to previous variables
    Sol.uPrev = Sol.u ;
    Sol.vPrev = Sol.v ;
    Sol.vDotPrev = Sol.vDot;  
    
    % Post-process the results
    % Write the integrated values of Metric and the forces in *.oisd file.
%     if (mod(timeStep,solver.intgOutFreq)==0)
%        fprintf(fileId1,'timeStep %d\n',timeStep);
%        fprintf(fileId1,'Metric\n');
%        fprintf(fileId1,'%e\n',Length);
%        fprintf(fileId1,'Traction\n');
%        fprintf(fileId1,'%e %e\n',Force(1),Force(2));
%        fprintf(fileId1,'Rigid_body_Displacement\n');
%        fprintf(fileId1,'%e %e\n',Sol.dispS(1),Sol.dispS(2));
%     end
%     
    % Write the results in Tecplot format (Change it to any other format or plot it in MATLAB)
    if (mod(timeStep,solver.outFreq)==0)
        filename = sprintf('./Transient_LinElast/%s.%d',problemString,timeStep);
        StoreVTK_Data(filename,problemType,ndof,nElem,nen,Sol.crd,cnn,Sol.u,Sol.v);
    end
    
end
fclose(fileId1);




