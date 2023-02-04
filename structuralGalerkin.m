function [LHS, RHS] = structuralGalerkin(Sol, solid, pmc, solver, ...
                        iif, jjf, fi, xxf, yyf, ...
                        ux, uy, gW, N, Nx, Ny, ...
                        nElem, nQuad, nen, ndof)


% accFac = (1.0-pmc.alphaM)/(1.0-pmc.alpha)*(1.0/(pmc.beta*solver.dt^2));     
% dispFac = pmc.beta*solver.dt/pmc.gamma;

accFac = pmc.alphaM/(pmc.alpha * (pmc.gamma1*pmc.gamma2));

lambda = solid.lambda;
mu = solid.mu;

                    
sA  = zeros(nen^2*nElem,nQuad);
sA1 = zeros(nen^2*nElem,nQuad); 
sA2 = zeros(nen^2*nElem,nQuad);
sA3 = zeros(nen^2*nElem,nQuad);

sA4 = zeros(nen^2*nElem,nQuad);
sA5 = zeros(nen^2*nElem,nQuad);
sA6 = zeros(nen^2*nElem,nQuad);

sA7 = zeros(nen*nElem,nQuad);
sA8 = zeros(nen*nElem,nQuad);
sA9 = zeros(nen*nElem,nQuad);
sA10 = zeros(nen*nElem,nQuad);

for p = 1:nQuad  
    J = [xxf*[Nx(:,p)], xxf*[Ny(:,p)],...
         yyf*[Nx(:,p)], yyf*[Ny(:,p)]];
    if size(J,2)==1
        J = J';
    end
    volume =( J(:,1).*J(:,4) -J(:,2).*J(:,3) );
           
    volume = abs(volume);
    
    negJacobian = find(volume<0);
    if ~isempty(negJacobian)
       disp('Mesh deformed, Negative Jacobian');
       exit
    end

    DNDx = ((J(:,4))*Nx(:,p)'+ (-J(:,3))*Ny(:,p)')./repmat(volume,1,nen);
    DNDy = ((-J(:,2))*Nx(:,p)'+(J(:,1))*Ny(:,p)')./repmat(volume,1,nen);   
            
    index = 0;
    index1 = 0;
    for i = 1:nen
        for j = 1:nen
             % Galerkin inertia term
            Mij = gW(p)*(N(p,i)*N(p,j));
            Mij = Mij.*solid.dens.*volume ;
            sA(index+1:index+nElem,p) = Mij;
            
            Aij_1 = gW(p)*DNDx(:,i).*DNDx(:,j);
            Aij_2 = gW(p)*DNDx(:,i).*DNDy(:,j);
            Aij_3 = gW(p)*DNDy(:,i).*DNDx(:,j);
            Aij_4 = gW(p)*DNDy(:,i).*DNDy(:,j);
            
            Aij_5 = gW(p)*N(p,i).*N(p,j);
            
            sA1(index+1:index+nElem,p) = ((lambda+2*mu)*Aij_1 + mu*Aij_4).*volume;
            sA2(index+1:index+nElem,p) = (lambda*Aij_2 + mu*Aij_3).*volume;
            sA3(index+1:index+nElem,p) = (lambda*Aij_3 + mu*Aij_2).*volume;
            sA4(index+1:index+nElem,p) = ((lambda+2*mu)*Aij_4 + mu*Aij_1).*volume;
            
%             sA5(index+1:index+nElem,p) = (Aij_1 + Aij_4).*volume;
            
            sA6(index+1:index+nElem,p) = Aij_5.*volume;
            
            index = index + nElem;
        end
        
%         F11_i = gW(p).*DNDx(:,i).*ux(:,i).*volume;
%         F12_i = gW(p).*DNDy(:,i).*ux(:,i).*volume;
%         F21_i = gW(p).*DNDx(:,i).*uy(:,i).*volume;
%         F22_i = gW(p).*DNDy(:,i).*uy(:,i).*volume;
        
%         sA7(index1+1:index1+nElem,p) = F11_i;
%         sA8(index1+1:index1+nElem,p) = F12_i;
%         sA9(index1+1:index1+nElem,p) = F21_i;
%         sA10(index1+1:index1+nElem,p) = F22_i;
%         index1 = index1+nElem;
    end
end
% Summation of all quadrature data
sA  = sum(sA, 2);
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
sA4 = sum(sA4,2);
% sA5 = sum(sA5,2);
sA6 = sum(sA6,2);
% sA7 = sum(sA7,2);
% sA8 = sum(sA8,2);
% sA9 = sum(sA9,2);
% sA10 = sum(sA10,2);
 
% Assemble the matrix
ZeroF = sparse(ndof,ndof);
        
Mf  = sparse(iif,jjf,sA, ndof,ndof);
A1 = sparse(iif,jjf,sA1,ndof,ndof);
A2 = sparse(iif,jjf,sA2,ndof,ndof);
A3 = sparse(iif,jjf,sA3,ndof,ndof);
A4 = sparse(iif,jjf,sA4,ndof,ndof);
% A5 = sparse(iif,jjf,sA5,ndof,ndof);
A6 = sparse(iif,jjf,sA6,ndof,ndof);

% F11 = 1 + sparse(fi,ones(size(fi)),sA7,ndof,ndof)*[ones(ndof,1)];
% F12 = sparse(fi,ones(size(fi)),sA8,ndof,ndof)*[ones(ndof,1)];
% F21 = sparse(fi,ones(size(fi)),sA9,ndof,ndof)*[ones(ndof,1)];
% F22 = 1 + sparse(fi,ones(size(fi)),sA10,ndof,ndof)*[ones(ndof,1)];

% E11 = 0.5*(F11.*F11 + F21.*F21 - 1.0);
% E12 = 0.5*(F11.*F12 + F21.*F22);
% E21 = 0.5*(F12.*F11 + F22.*F21);
% E22 = 0.5*(F12.*F12 + F22.*F22 - 1.0);
% 

Mf = [Mf ZeroF; ZeroF Mf];

% Kf1 = [A1 A2 ; A3 A4];
% Kf2 = [A5 ZeroF; ZeroF A5];
% Kf3 = [A1 A3 ; A2 A4];
% Kf = solid.lambda*Kf1 + solid.mu*(Kf2 + Kf3);
Kf = [A1 A2; A3 A4];

Src = [A6.*solid.gravFrc(1) ZeroF; ZeroF A6.*solid.gravFrc(2)];
   
% Left-hand side matrix
LHS = accFac*Mf + Kf*solver.dt^2;
% LHS = Kf;
        
clear sA1 sA2 sA3 Aij_1 Aij_2 Aij_3 sA4 sA5  sA6
clear Aij_4 Aij_5 Aij_6
clear A1 A2 A3 A4 A5 Kf1 Kf2 Kf3

uAlpha = Sol.uAlpha(:,:,1) ;
uAlpha = [uAlpha(:)];
vDotAlpha = Sol.vDotAlpha(:,:,1) ;
vDotAlpha = [vDotAlpha(:)];

gravVec = [ones(2*ndof,1)];

% Right-hand side vector
RHS = -(Mf * vDotAlpha(:) + Kf*uAlpha(:) - solid.dens*Src * gravVec)*solver.dt^2 ;
% RHS = -Kf*uAlpha + solid.dens*Src*gravVec;


end

                                      
                                  