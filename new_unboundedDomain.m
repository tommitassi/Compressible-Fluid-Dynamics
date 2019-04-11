%% Method of Characteristics - Acoustics 1D

%% 2. Creation of the mesh and computation of the solution - Unbounded domain

function [ x, u ] = new_unboundedDomain( x_in, u_in, c )

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

k = 1;
m = 1;
while ( n+k <= n*(n+1)/2 )
    
    for j = m+1 : n
        
        [x(:, n+k), u(:, n+k)] = innerProblem( x_in(:, m), x_in(:, j), ...
            u_in(:, m), u_in(:, j), c);
        k = k + 1;
        
    end
    
    m = m + 1;
    
end

end