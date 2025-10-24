% Poisson equation over a rectangular domain
% Using linear triangular finite elements
% Dirichlet boundary condition will be applied later

    % --- Create rectangular mesh ---
    a = 0; b = 2;      % x-range
    c = 0; d = 3;      % y-range
    nx_div = 50; ny_div = 50;   % number of divisions in x and y
    hx = (b - a) / nx_div;
    hy = (d - c) / ny_div;
    x = a:hx:b;
    y = c:hy:d;
    [X, Y] = meshgrid(x, y);
    nodes = [X(:), Y(:)];
    nNodes = size(nodes, 1);

    q0 = 500;        % source term
    h = 50;        % therman conductivity

    % --- Create triangular elements (2 per rectangle) ---
    nx = length(x);
    ny = length(y);
    elements = [];
    for j = 1:ny-1
        for i = 1:nx-1
            n1 = (j-1)*nx + i;       % bottom-left
            n2 = n1 + 1;             % bottom-right
            n3 = n1 + nx + 1;        % top-right
            n4 = n1 + nx;            % top-left
            elements = [elements;
                        n1 n2 n3;    % first triangle
                        n3 n4 n1];   % second triangle
        end
    end
    nElements = size(elements, 1);

    % --- Initialize global stiffness matrix and load vector ---
    K = zeros(nNodes, nNodes);
    F = zeros(nNodes, 1);

    % --- Loop over elements: compute local and assemble ---
    for e = 1:nElements
        conn = elements(e, :);           % 3 node indices
        coords = nodes(conn, :);         % coordinates of the 3 nodes
        [Ke, Fe] = local_stiffness(coords, q0, h);  % local matrices

        % Check if local stiffness matrix is symmetric
        if ~checkSymmetry(Ke)
            fprintf('Local stiffness matrix of element %d is NOT symmetric.\n', e);
        end

        % Assembly into global system
        for a = 1:3
            A = conn(a);
            F(A) = F(A) + Fe(a);
            for b = 1:3
                B = conn(b);
                K(A, B) = K(A, B) + Ke(a, b);
            end
        end
    end

    % --- Check global stiffness matrix symmetry ---
    if checkSymmetry(K)

        nx = nx_div + 1;  % number of nodes along x
        ny = ny_div + 1;  % number of nodes along y
        
        % left boundary: first column of nodes
        left = 1:nx:nNodes;
        
        % right boundary: last column of nodes
        right = nx:nx:nNodes;
        
        % bottom boundary: first row of nodes
        bottom = 1:nx;
        
        % top boundary: last row of nodes
        top = nNodes-nx+1 : nNodes;
    
        dirichletNodes = unique([left, right, bottom, top]);
    
    
        u_D = 100*ones(length(dirichletNodes), 1);
    
        % Initialize solution vector
        u = zeros(nNodes,1);
        
        % Apply Dirichlet conditions
        for k = 1:length(dirichletNodes)
            i = dirichletNodes(k);
            u(i) = u_D(k);
            K(i,:) = 0;           % zero out row
            K(i,i) = 1;           % put 1 on diagonal
            F(i) = u(i);          % enforce value
        end
        if ~checkSymmetry(K)
            % Solve
            u = K \ F;        
            % --- Output---        
            disp('Max U');
            disp(max(u))
        
            disp('Min U');
            disp(min(u))
        
            %figure
            %trisurf(elements, nodes(:,1), nodes(:,2), u, 'EdgeColor','k'); 
            %xlabel('x'); 
            %ylabel('y'); 
            %zlabel('u'); 
            %title('FEM solution with Dirichlet boundary');
        else
            disp('Global stiffness matrix is symmetric.');
        end

    else
        disp('Global stiffness matrix is NOT symmetric.');
    end

% -------------------------------------------------------------------------
function [Ke, Fe] = local_stiffness(N, q0, h)
    % N: 3x2 matrix of node coordinates of one triangle
    i = N(1, :);
    j = N(2, :);
    k = N(3, :);

    bi = j(2) - k(2);
    bj = k(2) - i(2);
    bk = i(2) - j(2);
    ci = k(1) - j(1);
    cj = i(1) - k(1);
    ck = j(1) - i(1);

    A = 0.5 * abs(det([1, i(1), i(2);
                       1, j(1), j(2);
                       1, k(1), k(2)]));

    K = [bi*bi + ci*ci,  bi*bj + ci*cj,  bi*bk + ci*ck;
         bj*bi + cj*ci,  bj*bj + cj*cj,  bj*bk + cj*ck;
         bk*bi + ck*ci,  bk*bj + ck*cj,  bk*bk + ck*ck];

    Ke = (h / (4*A)) * K;
    Fe = ((q0 * A) / 3) * ones(3, 1);
end

function isSymmetric = checkSymmetry(A)
    tol = 1e-12; % tolerance for floating-point comparisons
    isSymmetric = norm(A - A', 'fro') < tol;
end

