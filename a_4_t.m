% Clear the workspace
clear
clc

% Set the size of the matrix
Nx = 400;
Ny = Nx;

% Initialize the matrix A
A = zeros(Nx, Ny);

% Calculate the square root of Nx
s = sqrt(Nx);

% Fill the matrix A according to the problem's pattern
for i = 1:Nx
    for j = 1:Ny
        if i == j
            A(i, j) = -4;
        end
    end
end

for i = 1:Nx - s
    for j = s + 1:Ny
        if i == j - s
            A(i, j) = 1;
        end
    end
end

for i = s + 1:Nx
    for j = 1:Ny - s
        if i == j + s
            A(i, j) = 1;
        end
    end
end

for i = 2:Nx
    for j = 1:Ny - 1
        if i == j + 1 && mod(j, s) ~= 0
            A(i, j) = 1;
        end
    end
end

for i = 1:Nx - 1
    for j = 2:Ny
        if i == j - 1 && mod(i, s) ~= 0
            A(i, j) = 1;
        end
    end
end

% Initialize the right-hand side vector b
b = zeros(Nx, 1);
for i = 1:s
    b(i) = -1;
end

% Solve the linear system A\b using Gaussian elimination
E = A\b;

% Initialize variables for Jacobi iteration
x0 = zeros(Nx, 1);
max_iterations = 10000;
tolerance = 1e-6;

% Perform Jacobi iteration
x1 = x0;
for k = 1:max_iterations
    x_new = zeros(Nx, 1);
    for i = 1:Nx
        sum = 0;
        for j = 1:Nx
            if i ~= j
                sum = sum + A(i, j) * x1(j);
            end
        end
        x_new(i) = (1 / A(i, i)) * (b(i) - sum);
    end
    
    % Check for convergence
    if norm(x_new - x1) < tolerance
        break;
    end
    
    x1 = x_new;
end

fprintf("%d iterations in Jacobi Method\n", k);

% Initialize variables for Gauss-Seidel iteration
x0 = zeros(Nx, 1);
max_iterations = 10000;
tolerance = 1e-6;

% Perform Gauss-Seidel iteration
x2 = x0;
for k = 1:max_iterations
    for i = 1:Nx
        sum = 0;
        for j = 1:Nx
            if i ~= j
                sum = sum + A(i, j) * x2(j);
            end
        end
        x2(i) = (1 / A(i, i)) * (b(i) - sum);
    end
    
    % Check for convergence
    if norm(A * x2 - b) < tolerance
        break;
    end
end

fprintf("%d iterations in Gauss-Seidel Method\n", k);

% Initialize variables for SOR iteration
x0 = zeros(Nx, 1);
max_iterations = 10000;
tolerance = 1e-6;
w = 1.25; % Relaxation parameter 

% Perform SOR iteration
x3 = x0;
for k = 1:max_iterations
    for i = 1:Nx
        sum1 = 0;
        sum2 = 0;
        for j = 1:i - 1
            sum1 = sum1 + A(i, j) * x3(j);
        end
        for j = i + 1:Nx
            sum2 = sum2 + A(i, j) * x3(j);
        end
        x3(i) = (1 - w) * x3(i) + (w / A(i, i)) * (b(i) - sum1 - sum2);
    end
    
    % Check for convergence
    if norm(A * x3 - b) < tolerance
        break;
    end
end

fprintf("%d iterations in SOR Method\n", k);

% Reshape the solutions and prepare for plotting
E = reshape(E, [s, s]);
Z = zeros(s + 2, s + 2);
E = E';

for i = 2:s + 1
    for j = 2:s + 1
        Z(i, j) = E(i - 1, j - 1);
    end
end

Z(1, :) = 1;

x1 = reshape(x1, [s, s]);
Z1 = zeros(s + 2, s + 2);
x1 = x1';

for i = 2:s + 1
    for j = 2:s + 1
        Z1(i, j) = x1(i - 1, j - 1);
    end
end

Z1(1, :) = 1;

x2 = reshape(x2, [s, s]);
Z2 = zeros(s + 2, s + 2);
x2 = x2';

for i = 2:s + 1
    for j = 2:s + 1
        Z2(i, j) = x2(i - 1, j - 1);
    end
end

Z2(1, :) = 1;

x3 = reshape(x3, [s, s]);
Z3 = zeros(s + 2, s + 2);
x3 = x3';

for i = 2:s + 1
    for j = 2:s + 1
        Z3(i, j) = x3(i - 1, j - 1);
    end
end

Z3(1, :) = 1;

% Define the x and y coordinates for plotting
x = linspace(0, 1, s + 2);
y = linspace(1, 0, s + 2);
[X, Y] = meshgrid(x, y);


% Plot the results
figure;
subplot(1, 3, 1);
contourf(X, Y, Z1, 12, '--');
colorbar;
title('Jacobi Iteration');

subplot(1, 3, 2);
contourf(X, Y, Z2, 12, '--');
colorbar;
title('Gauss-Seidel Iteration');

subplot(1, 3, 3);
contourf(X, Y, Z3, 12, '--');
colorbar;
title('SOR Iteration');

sgtitle('Steady-State Solutions Comparison');
