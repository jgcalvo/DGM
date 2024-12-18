function [A,b,u,err] = DGM(mesh,rhs,alpha,uex,varargin)
% DGM computes the stiffness matrix (with boundary data) AD, the right-hand
% side bD, the exact solution u, and the L2-error err, for Poisson equation
% in 2D, given by -div(\rho \grad u) = rhs,
% discretized with the symmetric interior penalty discontinuous
% Galerkin element method on a polygonal mesh.
%
% SYNOPSIS: DGM(mesh,rhs,alpha,uex,rho)
%
% INPUT:  mesh:	 struct with the following fields:
%                mesh.nodes2coord - mesh nodes
%                mesh.elems2nodes - mesh connectivity
%                mesh.edges       - edge connectivity
%                mesh.edges2elems - edge2elem matrix
%                mesh.edgeSign    - orientation of edges
%                mesh.NE          - number of elements
%                mesh.Ne          - number of edges
%                mesh.NeBnd       - number of boundary edges
%                mesh.NeInt       - number of interior edges
%         rhs:   function handle for the right-hand side
%         alpha: stability constant
%         uex:   function handle for the exact solution (to cumpute err)
%         rho:   function handle for the coefficient
%
% EXAMPLE:
%         DGM(mesh,rhs,alpha,uex) uses rho = 1

% AUTHOR: Juan G. Calvo and collaborators, 2023

% read fields of mesh
nodes2coord = mesh.nodes2coord;
elems2nodes = mesh.elems2nodes;
edges2elems = mesh.edges2elems;
edges       = mesh.edges;
edgeSign    = mesh.edgeSign;
NE    = mesh.NElems;
Ne    = mesh.Nedges;
NeBnd = mesh.NedgesBnd;
NeInt = mesh.NedgesInt;

% coefficient $\rho$
if(size(varargin,2)==1)  % read \rho (if given)
    rhoFn    = varargin{1};
    % evaluation of \rho at centroids of triangles
    centroid = (nodes2coord(elems2nodes(:,1),:)+nodes2coord(elems2nodes(:,2),:)+nodes2coord(elems2nodes(:,3),:))/3;
    rhoElem  = rhoFn(centroid(:,1),centroid(:,2));
    % evaluation of \rho at midpoints of edges
    midptsEdges = (nodes2coord(edges(:,1),:)+nodes2coord(edges(:,2),:))/2;
    rhoEdge = rhoFn(midptsEdges(:,1),midptsEdges(:,2));
else                     % otherwise \rho = 1
    rhoElem = ones(size(elems2nodes,1),1);
    rhoEdge = ones(size(edges,1),1);
end

% contributions of \int_K \rho \nabla \phi_i \cdot \nabla \phi_j
glbInd = 1:3*NE; glbInd = reshape(glbInd,3,NE)';
v1 = nodes2coord(elems2nodes(:,1),:);
v2 = nodes2coord(elems2nodes(:,2),:);
v3 = nodes2coord(elems2nodes(:,3),:);
xx = v1-v2; yy = v2-v3; zz = v3-v1;
area = 0.5*(-xx(:,1).*zz(:,2) + xx(:,2).*zz(:,1));
a11 = sum(yy.*yy,2)./(4*area).*rhoElem;
a12 = sum(yy.*zz,2)./(4*area).*rhoElem;
a13 = sum(yy.*xx,2)./(4*area).*rhoElem;
a22 = sum(zz.*zz,2)./(4*area).*rhoElem;
a23 = sum(zz.*xx,2)./(4*area).*rhoElem;
a33 = sum(xx.*xx,2)./(4*area).*rhoElem;
ii = [glbInd(:,1); glbInd(:,1); glbInd(:,2); glbInd(:,1); glbInd(:,3); glbInd(:,2); glbInd(:,2); glbInd(:,3); glbInd(:,3)];
jj = [glbInd(:,1); glbInd(:,2); glbInd(:,1); glbInd(:,3); glbInd(:,1); glbInd(:,2); glbInd(:,3); glbInd(:,2); glbInd(:,3)];
ss = [a11; a12; a12; a13; a13; a22; a23; a23; a33];
A1 = sparse(ii,jj,ss,3*NE,3*NE);
clear a11 a12 a13 a22 a23 a33 ii jj ss v1 v2 v3 xx yy zz

% contributions of edge integrals
% allocate vectors
sizeMat = 9*NeBnd + 36*NeInt; pt = 0;      % matrices: 3x3 edge on boundary, 6x6 interior edge
ii = nan(sizeMat,1); ssB = nan(sizeMat,1); % allocate space
jj = nan(sizeMat,1); ssC = nan(sizeMat,1); % for matrices
for edgeID = 1:Ne                    % loop over all edges
    elem1 = edges2elems(edgeID,1);   % read element K+
    elem2 = edges2elems(edgeID,2);   % read element K-
    if(elem1-elem2~=0)               % interior edge
        vP = elems2nodes(elem1,myismember(elems2nodes(elem1,:),edges(edgeID,:)));
        kN = elems2nodes(elem2,myismember(elems2nodes(elem2,:),edges(edgeID,:)));
        v1  = nodes2coord(edges(edgeID,1),:);
        v2  = nodes2coord(edges(edgeID,2),:);
        v3p = nodes2coord(vP,:);
        v3n = nodes2coord(kN,:);
        % check orientation
        if(edgeSign(edgeID)<0)       % fix orientation
            v3p = nodes2coord(kN,:);
            v3n = nodes2coord(vP,:);
            tem = elem2;
            elem2 = elem1;
            elem1 = tem;
        end
        xx  = v2-v1; yy = v3p-v2; zz = v1-v3p; aa = v3n-v2; bb = v1-v3n;
        b11 = dot(yy,xx)/(4*mycross(yy,xx));
        b21 = dot(xx,zz)/(4*mycross(xx,zz));
        b31 = dot(xx,xx)/(4*mycross(yy,xx));
        B11 = [b11 b11 0; b21 b21 0; b31 b31 0];
        b11 = dot(xx,aa)/(4*mycross(xx,aa));
        b21 = dot(bb,xx)/(4*mycross(bb,xx));
        b31 = dot(xx,xx)/(4*mycross(xx,aa));
        B22 = [b11 b11 0; b21 b21 0; b31 b31 0];
        % bilinear form b1
        %M1 = [B11*rho(elem1) -B11*rho(elem2); -B22*rho(elem1) B22*rho(elem2)];
        M1 = rhoEdge(edgeID)*[B11 -B11; -B22 B22];
        % bilinear form b2
        B11 = [1/3 1/6 0; 1/6 1/3 0; 0 0 0];
        M2  = [B11 -B11; -B11 B11];
        % position: find global dofs
        ind1P = elems2nodes(elem1,:)==edges(edgeID,1);
        ind2P = elems2nodes(elem1,:)==edges(edgeID,2);
        ind3P = find(~(ind1P+ind2P));
        ind1P = find(ind1P);
        ind2P = find(ind2P);
        loc1  = glbInd(elem1,:);
        loc1  = loc1([ind1P ind2P ind3P]);
        ind1N = elems2nodes(elem2,:)==edges(edgeID,1);
        ind2N = elems2nodes(elem2,:)==edges(edgeID,2);
        ind3N = find(~(ind1N+ind2N));
        ind1N = find(ind1N);
        ind2N = find(ind2N);
        loc2  = glbInd(elem2,:);
        loc2  = loc2([ind1N ind2N ind3N]);
        [yy,xx] = meshgrid([loc1 loc2]);
        % store data
        adv = 36;
        ii(pt+1:pt+adv)  = xx(:);
        jj(pt+1:pt+adv)  = yy(:);
        ssB(pt+1:pt+adv) = M1(:);
        ssC(pt+1:pt+adv) = M2(:);
        pt = pt + adv;
    else                             % boundary edge
        vP  = setdiff(elems2nodes(elem1,:),edges(edgeID,:));
        % set orientation
        if(edgeSign(edgeID)==1)
            v1  = nodes2coord(edges(edgeID,1),:);
            v2  = nodes2coord(edges(edgeID,2),:);
        else
            v1  = nodes2coord(edges(edgeID,2),:);
            v2  = nodes2coord(edges(edgeID,1),:);
        end
        v3p = nodes2coord(vP,:);
        xx  = v2-v1; yy = v3p-v2; zz = v1-v3p;
        b1 = dot(yy,xx)/(2*mycross(yy,xx));
        b2 = dot(xx,zz)/(2*mycross(xx,zz));
        b3 = dot(xx,xx)/(2*mycross(yy,xx));
        % bilinear form b1
        M1 = rhoEdge(edgeID)*[b1 b1 0; b2 b2 0; b3 b3 0];
        % bilinear form b2
        M2 = [1/3 1/6 0; 1/6 1/3 0; 0 0 0];
        % position: find global dofs
        ind1P = elems2nodes(elem1,:)==edges(edgeID,1);
        ind2P = elems2nodes(elem1,:)==edges(edgeID,2);
        ind3P = find(~(ind1P+ind2P));
        ind1P = find(ind1P);
        ind2P = find(ind2P);
        loc1  = glbInd(elem1,:);
        loc1  = loc1([ind1P ind2P ind3P]);
        [yy,xx] = meshgrid(loc1);
        % store data
        adv = 9;
        ii(pt+1:pt+adv) = xx(:);
        jj(pt+1:pt+adv) = yy(:);
        ssB(pt+1:pt+adv) = M1(:);
        ssC(pt+1:pt+adv) = M2(:);
        pt = pt + adv;
    end
end
M1 = sparse(ii,jj,ssB);
M2 = sparse(ii,jj,ssC);

% assemble matrix
A = A1-M1-M1'+alpha*M2;
clear A1 M1 M2

% rhs
% one point quadrature rule for rhs (centroid)
b = zeros(NE,3);
pxy = (nodes2coord(elems2nodes(:,1),:) ...
    + nodes2coord(elems2nodes(:,2),:) ...
    + nodes2coord(elems2nodes(:,3),:)) / 3;
for i = 1:3
    b(:,i) = b(:,i) + rhs(pxy)/3;
end
b = (b.*repmat(area,1,3))'; b = b(:);

% direct solver
u = A\b;

% compute L2 error
uElem = reshape(u,3,NE); uElem = uElem';
err = getL2error(nodes2coord,elems2nodes,uex,uElem,area);
end

%% required functions
function v = mycross(v1,v2)
v = v1(1)*v2(2)-v1(2)*v2(1);
end

function lia = myismember(a,b)
% find missing index of b in vector a
% Example: a = [4 6 7];  b = [6 4]; lia = 3
lia = 0;
add = 1;
while(add>0)
    lia = lia + 1;
    add = sum(b==a(lia));
end
end

function L2err = getL2error(node,elem,uex,uElem,area)
% compute L2err by adding L2 approximations of integral with one-point
% quadrature rule (evaluation at barycenter)
% approximated solution at barycenters
uh = sum(uElem,2)/3; 
% exact solution at barycenters
centroid = (node(elem(:,1),:)+node(elem(:,2),:)+node(elem(:,3),:))/3;
% compute quadrature rule
L2err = sqrt(sum(area.*(uex(centroid)-uh).^2));
end