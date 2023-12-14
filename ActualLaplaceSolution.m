






N = 10; % Number of terms in the series
L = 1;  % Length of the domain in the x-direction
H = 1;  % Height of the domain in the y-direction
g = @(y) 4 .* y .* (H - y); % Corrected: Boundary Condition as a function

laplacianFourierSolution(N, L, H, g); % Pass g as an argument

function laplacianFourierSolution(N, L, H, g)
    % N: Number of terms in the series
    % L: Length of the domain in the x-direction
    % H: Height of the domain in the y-direction

    % Create a grid for plotting
    [X, Y] = meshgrid(linspace(0, L, 100), linspace(0, H, 100));

    % Initialize the solution u(x, y) to zero
    U = zeros(size(X));
    
    % Compute the first N terms of the series
    for n = 1:N
        % Define the function for the integral
        fun = @(y) g(y) .* sin(n * pi * y / H);
        
        % Compute the Fourier coefficient Bn
        Bn = 2 * integral(fun, 0, H) / (H * sinh(n * pi * L / H));
        
        % Update the solution U
        U = U + Bn * sin(n * pi * Y / H) .* sinh(n * pi * X / H);
    end
    
    % Define the query point
    xQuery = 0.412; yQuery = 0.563; % Test point

    % Interpolate the value of U at the query point
    uQuery = interp2(X, Y, U, xQuery, yQuery);


    % Plot the solution
    figure('WindowStyle', 'docked');
    surf(X, Y, U);
    title(['Solution of the 2D Laplacian, First ' num2str(N) ' Fourier Terms']);
    xlabel('x');
    ylabel('y');
    zlabel('u(x, y)');
    colorbar;

    % Display the interpolated value at the query point
    disp(['u(', num2str(xQuery), ', ', num2str(yQuery), ') = ', num2str(uQuery)]);
    disp(' ')
end
