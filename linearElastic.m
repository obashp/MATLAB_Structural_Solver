%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %  
%              2D Hookean Solid Finite-element Solver                     %
%                                                                         % 
%      \rho u_tt = div(C sym(\grad u)) + \rho g  in \Omega,               %
%                                                                         %
%       Dirichlet boundary condition        u = g_D  on \Gamma_D,         %
%       Neumann boundary condition du/dn - np = g_N  on \Gamma_N.         %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sol, HEnormIndicator] = linearElastic(solver, solid, pmc, Sol, cnn, crd, ...
                                               elemType, ndof, nen, nElem)

% Quadrature rules for elements
if strcmp(elemType,'3Tri')
    gP = ...
   [1/3,  1/3
    4/3,  1/3
    1/3,  4/3] ;
    gW = ...
   [2/3,  2/3,  2/3] ;
 
    N(:,1) = 0.5.*(2.-gP(:,1)-gP(:,2)) ;
    N(:,2) = 0.5.*(gP(:,1)) ;
    N(:,3) = 0.5.*(gP(:,2)) ;
    
    Nx(:,1) = -0.5.*ones(3,1) ;
    Nx(:,2) =  0.5.*ones(3,1) ;
    Nx(:,3) =  zeros(3,1) ; 
    Ny(:,1) = -0.5.*ones(3,1) ;
    Ny(:,2) =  zeros(3,1) ;
    Ny(:,3) =  0.5.*ones(3,1) ;    
elseif strcmp(elemType,'4Quad')
    gP = ...
   [-5.7735026918962584E-01, -5.7735026918962584E-01
     5.7735026918962584E-01, -5.7735026918962584E-01
     5.7735026918962584E-01,  5.7735026918962584E-01
    -5.7735026918962584E-01,  5.7735026918962584E-01] ;
    gW = [1, 1, 1, 1 ] ;
    
    N(:,1) = 0.25.*(1-gP(:,1)).*(1-gP(:,2)) ;
    N(:,2) = 0.25.*(1+gP(:,1)).*(1-gP(:,2)) ;
    N(:,3) = 0.25.*(1+gP(:,1)).*(1+gP(:,2)) ;
    N(:,4) = 0.25.*(1-gP(:,1)).*(1+gP(:,2)) ;
    
    Nx(:,1) = -0.25.*(1-gP(:,2)) ;
    Nx(:,2) =  0.25.*(1-gP(:,2)) ;
    Nx(:,3) =  0.25.*(1+gP(:,2)) ;
    Nx(:,4) = -0.25.*(1+gP(:,2)) ;
    Ny(:,1) = -0.25.*(1-gP(:,1)) ;
    Ny(:,2) = -0.25.*(1+gP(:,1)) ;
    Ny(:,3) =  0.25.*(1+gP(:,1)) ;
    Ny(:,4) =  0.25.*(1-gP(:,1)) ;
end

Nx = Nx' ;
Ny = Ny' ;
nQuad = length(gW) ;
 
iif = zeros(nen^2*nElem,1); 
jjf = zeros(nen^2*nElem,1);
index = 0;
for i = 1:nen
   for j = 1:nen
      iif(index+1:index+nElem) = double(cnn(:,i)); 
      jjf(index+1:index+nElem) = double(cnn(:,j));  
      index = index + nElem;
   end
end

index = 0;
for i = 1:nen
    fi(index+1:index+nElem) = double(cnn(:,i));
    index = index+nElem;
end
        
% Satisfy boundary conditions
Sol.u(solid.DirichletU,1,1) = solid.DirichletUval ;
Sol.u(solid.DirichletV,2,1) = solid.DirichletVval ;

% Sol.uAlpha = Sol.u;
% Sol.vAlpha = Sol.v;
        
% Navier-Stokes equations
xxf = zeros(size(cnn));
yyf = zeros(size(cnn));
ux = zeros(size(cnn));
uy = zeros(size(cnn));

for i=1:nen
   xxf(:,i) = crd(cnn(:,i),1);
   yyf(:,i) = crd(cnn(:,i),2);
   ux(:,i) = Sol.uAlpha(cnn(:,i),1,1);
   uy(:,i) = Sol.uAlpha(cnn(:,i),2,1);
end
        
% Form element matrix and assemble Galerkin terms
[LHS,RHS] = structuralGalerkin(Sol, solid, pmc, solver, iif, jjf, fi, xxf, yyf, ...
                        ux, uy, gW, N, Nx, Ny, ...
                        nElem, nQuad, nen, ndof);
                                      
% Solve the linear system

% Select the unknown nodal values
freeNodesUx = unique([solid.DirichletU]);
freeNodesUx = setdiff(1:size(crd,1),[freeNodesUx]);
freeNodesUy = unique([solid.DirichletV]);
freeNodesUy = setdiff(1:size(crd,1),[freeNodesUy]);

freeNodes = [freeNodesUx';freeNodesUy' + size(crd,1)];
        
result = Sol.uAlpha(:,:,1);
result = result(:);
% resultDot = Sol.vAlpha(:,:,1);
% resultDot = resultDot(:);
% resultDDot = Sol.vDotAlpha(:,:,1);
% resultDDot = resultDDot(:);

Increment = LHS(freeNodes,freeNodes)\RHS(freeNodes);
        
% Update the increments

% Alpha1 = (1.0-pmc.alphaM)/(1.0-pmc.alpha)*(1.0/(pmc.beta*solver.dt^2));
% Alpha2 = (pmc.alphaM -1.0)/(pmc.beta*solver.dt);
% Alpha3 = Alpha2*solver.dt*(0.5 - pmc.beta/(1.0-pmc.alphaM));

Alpha1 = pmc.alphaM/(pmc.alpha*pmc.gamma1*pmc.gamma2*solver.dt^2);
Alpha2 = -pmc.alphaM/(pmc.gamma1*pmc.gamma2*solver.dt);
Alpha3 = (1.0 - pmc.alphaM/pmc.gamma2);



result(freeNodes) = result(freeNodes) + Increment;
% resultDot(freeNodes) = resultDot(freeNodes) + (pmc.gamma/(pmc.beta*solver.dt))*Increment ;
% resultDDot(freeNodes) = resultDDot(freeNodes) ...
%        + (1.0-pmc.alphaM)/(1.0-pmc.alpha)*(1.0/(pmc.beta*solver.dt^2))*Increment;

% Sol.u(:,:,1) = reshape(result(1:2*ndof),[],2);
Sol.uAlpha(:,:,1) = reshape(result(1:2*ndof),[],2);
% Sol.vAlpha(:,:,1) = reshape(resultDot(1:2*ndof),[],2);
% Sol.vDotAlpha(:,:,1) = reshape(resultDDot(1:2*ndof),[],2);
% Sol.vDotAlpha =  Alpha1*(Sol.uAlpha-Sol.uPrev) + Alpha2*Sol.vPrev + Alpha3*Sol.vDotPrev;
Sol.vDotAlpha = Alpha1*(Sol.uAlpha - Sol.uPrev) + Alpha2*Sol.vPrev + Alpha3*Sol.vDotPrev;
Sol.vAlpha = (1 - pmc.alpha/pmc.gamma1)*Sol.vPrev + (1.0/(pmc.gamma1*solver.dt))*(Sol.uAlpha - Sol.uPrev);



% 
% ((pmc.alphaM - 1 + pmc.beta)/pmc.beta)*Sol.vDotPrev ...
%     + (1.0/(pmc.beta*solver.dt^2))...
%     * (((1.0-pmc.alphaM)/(1.0-pmc.alpha))*(Sol.uAlpha - Sol.uPrev) ...
%     - (1.0-pmc.alpha)*solver.dt*Sol.vPrev);
% 
% Sol.vAlpha = Sol.vPrev ...
%     + ((1.0-pmc.alpha)/(1.0-pmc.alphaM)*solver.dt)...
%     * (pmc.gamma*Sol.vDotAlpha + (1.0-pmc.alphaM-pmc.gamma)*Sol.vDotPrev);



HEnormIndicator =  norm(Increment)/norm(result(freeNodes)) ;
fprintf('HE: %e , ', HEnormIndicator);


clear freeNodes1 freeNodes2 freeNodes3
clear result resultDot
        
end

