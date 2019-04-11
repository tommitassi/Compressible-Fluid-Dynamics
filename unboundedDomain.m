%% Method of Characteristics - Acoustics 1D

%% 2. Creation of the mesh and computation of the solution - Unbounded domain

function [ x, u ] = unboundedDomain( x_in, u_in, c )

% unboundedDomain creates the mesh of the unbounded problem
% from the points where the initial conditions are imposed
% it creates the mesh in their domain of influence and evaluetes 
% the solution in those points

% INPUT
%   x_in    [2 x n]         Coordinates of the points where the initial
%                           condition are imposed x_in(:,i)=[x;t]
%   u_in    [2 x n]         Initial conditions u_in(:,i)=[v;w]
%   c       [1 x 1]         Velocity of propagation

% OUTPUT
%   x       [2 x n*(n+1)/2] Coordinates of the mesh points x(:,i)=[x;t]
%   u       [2 x n*(n+1)/2] Solution in mesh points u(:,i)=[v;w]

n = size(x_in, 2); % number of initial condition points
x = zeros(2, n*(n+1)/2);
u = zeros(2, n*(n+1)/2);

for i = 1 : n
    x(:, i) = x_in(:,i); % it copies the first step
    u(:, i) = u_in(:,i);
end

k = 1; % k is a parameter that allows to create the x vector in order
       % from left to right from step 0 to the last step
for i = n+1 : n*(n+1)/2
    
    [x(:, i), u(:,i)] = innerProblem( x(:, i-n+k-1), x(:, i-n+k), ...
        u(:, i-n+k-1), u(:, i-n+k), c);
    
    if i == (k+1)*n-(k+1)*k/2
        k = k + 1;
    end
    
end

end