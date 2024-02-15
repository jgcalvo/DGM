%% DGM
% parameters DGM
type = 'unstruc';                                  % 'struc' or 'unstruc'
uex = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));         % exact solution
rhs = @(x) 2*pi*pi*sin(pi*x(:,1)).*sin(pi*x(:,2)); % f
alpha = 10;                                        % stability constant
results = zeros(10,2); pt = 0;                     % store results
for N = 5%1:4                                        % loop over meshes
    % create triangulation
    mesh = load(['./meshes/tria_' type '_' num2str(N)]); mesh = mesh.mesh;
    % create stiffness matrix
    [AD,bD] = DGM(mesh,rhs,alpha);
    % direct solver
    u = AD\bD;
    % compute error (inf norm on nodes)
    uElem = reshape(u,3,mesh.NElems); uElem = uElem';
    uEx = nan(mesh.NElems,3);
    for j = 1:mesh.NElems
        uEx(j,:) = uex(mesh.nodes2coord(mesh.elems2nodes(j,:),:));
    end
    err = norm(uEx(:)-uElem(:),Inf);
    % store h and error
    results(pt+1,:) = [mesh.hmax err]; 
    pt = pt + 1;
end
% convergence plot
loglog(results(1:pt,1),results(1:pt,2),'.-')
% % plot solution
% plotSol(mesh,uElem), view(2), colorbar


function plotSol(mesh,solution)
nodes = mesh.nodes2coord;
elems = mesh.elems2nodes;
NE = size(elems,1);
for elID = 1:NE
    verts = nodes(elems(elID,:),:);
    trisurf([1 2 3], verts(:,1),verts(:,2), solution(elID,:)), hold on
end
end
