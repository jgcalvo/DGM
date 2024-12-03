function demo(opt,type)
% demo includes three examples for the solution of Poisson equation in 2D
% discretized with the symmetric interior penalty discontinuous Galerkin
% element method on a polygonal mesh
%
% SYNOPSIS: demo(opt,type)
%
% INPUT:  opt:	scalar with the following options:
%               1: convergence example as a function of h for a small
%               set of values of h (due to running times) for Poisson
%               equation
%               2: plot solution for a particular value of h
%               3: convergence example for rho = 1+xy
%         type: type of mesh; use 'struc' or 'unstruc' for structured
%                   or unstructured triangular meshes
%
% EXAMPLE:
%         demo(1,'struc'), demo(2,'unstruc'),

% AUTHOR: Juan G. Calvo and collaborators, 2023

close all

switch opt
    case {1,3} % convergence as a function of h for Poisson equation
        % exact solution
        uex = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));
        if(opt == 1) % case \rho = 1
            % rho coefficient for -nabla(rho \div u) = rhs
            rho = @(x,y) 1+0*x;
            % right hand side
            rhs = @(x) 2*pi*pi*sin(pi*x(:,1)).*sin(pi*x(:,2));
            % stability constant
            alpha = 10;
        elseif(opt == 3) % case \rho = 1+xy
            % rho coefficient for -nabla(rho \div u) = rhs
            rho = @(x,y) 1+x.*y;
            % right hand side
            rhs = @(x) -pi*x(:,2).*cos(pi*x(:,1)).*sin(pi*x(:,2)) ...
                       -pi*x(:,1).*sin(pi*x(:,1)).*cos(pi*x(:,2)) + ...
                2*pi*pi*sin(pi*x(:,1)).*sin(pi*x(:,2)).*rho(x(:,1),x(:,2));
            % stability constant
            alpha = 20;
        end
        results = zeros(10,2); % store results
        for N = 1:7            % loop meshes
            % load triangulation
            mesh = load(['./meshes/tria_' type '_' num2str(N)]); mesh = mesh.mesh;
            % assemble, solve and compute L2 error
            [~,~,~,err] = DGM(mesh,rhs,alpha,uex,rho);
            % store h and error
            results(N,:) = [mesh.hmax err];
        end
        %%
        % convergence plot
        loglog(results(:,1),results(:,2),'.-'), hold on
        loglog(results(:,1),results(:,1).^2,'.--'), hold off
        grid on
        xlabel('$h$','Interpreter','latex')
        ylabel('$\|u-u_h\|_{L^2}$','Interpreter','latex')
        set(gca,'FontSize',12)
        legend({'$error$','$h^2$'},'Interpreter','latex','Location','northwest')
    case 2 % plot solution for a particular choice of h
        % load triangulation
        uex = @(x) sin(pi*x(:,1)).*sin(pi*x(:,2));         % exact solution
        rhs = @(x) 2*pi*pi*sin(pi*x(:,1)).*sin(pi*x(:,2)); % f
        alpha = 10;                                        % stability constant
        mesh = load(['./meshes/tria_' type '_3']); mesh = mesh.mesh;
        % solve PDE
        [~,~,u,~] = DGM(mesh,rhs,alpha,uex);
        % reshape solution according to element matrix
        uElem = reshape(u,3,mesh.NElems); uElem = uElem';
        % exact solution at nodes
        uEx = nan(mesh.NElems,3);
        for j = 1:mesh.NElems
            uEx(j,:) = uex(mesh.nodes2coord(mesh.elems2nodes(j,:),:));
        end
        % plot solution and absolute error
        subplot(1,3,1)
        plotSol(mesh,uElem), shading interp
        title('Solution')
        subplot(1,3,2)
        plotSol(mesh,uElem), shading interp, view(2), colorbar
        title('Solution')
        subplot(1,3,3)
        plotSol(mesh,abs(uElem-uEx)), shading interp, view(2), colorbar
        title('Absolute Error')
end
end

function plotSol(mesh,solution)
% plotSol plots a discontinuous linear function for a given matrix with
% nodal values solution of the same size as the element matrix
nodes = mesh.nodes2coord; % read nodes
elems = mesh.elems2nodes; % read elements
NE = size(elems,1);       % number of elements
for elID = 1:NE           % plot linear function for each element (trisurf)
    verts = nodes(elems(elID,:),:);
    trisurf([1 2 3], verts(:,1),verts(:,2), solution(elID,:)), hold on
end
end
