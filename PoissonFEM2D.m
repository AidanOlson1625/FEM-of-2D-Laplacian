





FEM_2D_Poisson_Piecewise

function FEM_2D_Poisson_Piecewise
    % Square Domain
    xMin = 0; 
    yMin = 0;
    xMax = 1; 
    yMax = 1;
    
    numElemsX = 80; % Number of elements along x
    numElemsY = 80; % Number of elements along y

    % Boundary Conditions
    % top = @(x) 0; bottom = @(x) 0; left = @(y) 0; right = @(y) 0; % Dirichlet
    % top = @(x) 0; bottom = @(x) 0; left = @(y) 0; right = @(y) 4*y*(yMax-y); % Partially Dirichlet
    top = @(x) sin(4*pi*x + pi); bottom = @(x) sin(4*pi*x); left = @(y) sin(4*pi*y + pi); right = @(y) sin(4*pi*y); % Sinusoidal

    % Source function
    % f = @(x, y) exp(-((x-1/2).^2 +(y-1/2).^2)); % Gaussian
    % f = @(x, y) 800*((x-1/2).^2 - (y-1/2).^2); % Saddle
    f = @(x, y) 0; % Homogeneous

    % Generate mesh
    [nodeCoords, elementNodes] = generateMesh(xMin, xMax, yMin, yMax, numElemsX, numElemsY);

    % Solve FEM problem
    u = solveFEM(xMin, xMax, yMin, yMax, bottom, top, left, right, nodeCoords, elementNodes, f);

    % Query point
    xQuery = 0.412; yQuery = 0.563; % Test point for accuracy
    uQuery = evaluateSolution(xQuery, yQuery, nodeCoords, elementNodes, u);
    disp(['u(', num2str(xQuery), ', ', num2str(yQuery), ') = ', num2str(uQuery)]);

    % Construct and plot the piecewise linear solution
    plotPiecewiseSolution(numElemsX, numElemsY, nodeCoords, elementNodes, u);
    disp(['numElemsX = ' num2str(numElemsX) ' numElemsY = ' num2str(numElemsY)])
    disp(' ')
end

function [nodeCoords, elementNodes] = generateMesh(xMin, xMax, yMin, yMax, numElemsX, numElemsY)
    [X, Y] = meshgrid(linspace(xMin, xMax, numElemsX + 1), linspace(yMin, yMax, numElemsY + 1));
    nodeCoords = [X(:), Y(:)];
    
    elementNodes = delaunay(nodeCoords(:,1), nodeCoords(:,2));
end

function u = solveFEM(xMin, xMax, yMin, yMax, bottom, top, left, right, nodeCoords, elementNodes, f)
    numNodes = size(nodeCoords, 1);
    K = sparse(numNodes, numNodes); % Stiffness matrix
    F = zeros(numNodes, 1); % Load vector

    for e = 1:size(elementNodes, 1)
        nodeInd = elementNodes(e, :);
        x = nodeCoords(nodeInd, 1);
        y = nodeCoords(nodeInd, 2);

        % Element stiffness matrix and load vector
        [Ke, Fe] = elementStiffness(x, y, f);

        % Assemble
        K(nodeInd, nodeInd) = K(nodeInd, nodeInd) + Ke;
        F(nodeInd) = F(nodeInd) + Fe;
    end

    [K, F] = applyBC(K, F, nodeCoords, xMin, xMax, yMin, yMax, bottom, top, left, right);

    % Solve the system for constants
    u = K \ F;
end

function [K, F] = applyBC(K, F, nodeCoords, xMin, xMax, yMin, yMax, bottom, top, left, right)
    tol = 1e-10; % Tolerance
    for i = 1:size(nodeCoords, 1)
        x = nodeCoords(i,1);
        y = nodeCoords(i,2);
            if abs(x-xMin) < tol % Left boundary
                K(i,:) = 0;
                K(i,i) = 1;
                F(i) = left(y);
            elseif abs(x-xMax) < tol % Right boundary
                K(i,:) = 0;
                K(i,i) = 1;
                F(i) = right(y);
            elseif abs(y-yMin) < tol % Bottom boundary
                K(i,:) = 0;
                K(i,i) = 1;
                F(i) = bottom(x);
            elseif abs(y-yMax) < tol % Top boundary
                K(i,:) = 0;
                K(i,i) = 1;
                F(i) = top(x);
            end
    end
end

function [Ke, Fe] = elementStiffness(x, y, f)
    % Coordinates of the triangle vertices
    x1 = x(1); y1 = y(1);
    x2 = x(2); y2 = y(2);
    x3 = x(3); y3 = y(3);

    % Area of the triangle
    A = abs((x2-x1)*(y3-y1) - (x3-x1)*(y2-y1)) / 2;

    % Gradients of shape functions
    B = [y2-y3, y3-y1, y1-y2; x3-x2, x1-x3, x2-x1] / (2*A);

    % Element stiffness matrix
    Ke = A * (B' * B);

    % Element load vector - using centroid for approximation
    xc = mean(x); yc = mean(y);
    Fe = A * f(xc, yc) * [1; 1; 1] / 3;
end

function uQuery = evaluateSolution(xQuery, yQuery, nodeCoords, elementNodes, u)
    for e = 1:size(elementNodes, 1)
        nodeInd = elementNodes(e, :);
        x = nodeCoords(nodeInd, 1);
        y = nodeCoords(nodeInd, 2);

        % Barycentric coordinates for the query point
        detT = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));
        lambda1 = ((y(2)-y(3))*(xQuery-x(3)) + (x(3)-x(2))*(yQuery-y(3))) / detT;
        lambda2 = ((y(3)-y(1))*(xQuery-x(3)) + (x(1)-x(3))*(yQuery-y(3))) / detT;
        lambda3 = 1 - lambda1 - lambda2;

        % Check if the point is inside the triangle
        if lambda1 >= 0 && lambda1 <= 1 && lambda2 >= 0 && lambda2 <= 1 && lambda3 >= 0 && lambda3 <= 1
            uQuery = lambda1*u(nodeInd(1)) + lambda2*u(nodeInd(2)) + lambda3*u(nodeInd(3));
            return;
        end
    end
    error('Query point is outside the domain.');
end

function plotPiecewiseSolution(numElemsX, numElemsY, nodeCoords, elementNodes, u)
    figure('WindowStyle', 'docked')
    trisurf(elementNodes, nodeCoords(:,1), nodeCoords(:,2), u, 'FaceColor', 'interp');
    title('2D FEM Solution of the Poisson Equation');
    subtitle(['Number of Elements = ' num2str(2*numElemsX*numElemsY)]);
    xlabel('x'); ylabel('y'); zlabel('u(x, y)');
    colorbar;
end