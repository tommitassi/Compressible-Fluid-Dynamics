%% Method of Characteristics - Acoustics 1D

%% 3. Creation of the mesh relative to the left boundary

function [ x, u ] = leftBoundary_Still( x_in, u_in, ub, type, c )

% leftBoundary permits to complete the mesh using the left
% boundary condition and to compute the solution in those points

% INPUT
%   x_in    [2 x n]         Coordinates of the points where the initial
%                           condition are imposed x_in(:,i)=[x;t]
%   u_in    [2 x n]         Initial conditions u_in(:,i)=[v;w]
%   ub      [1 x 1]         Boundary condition
%   type    string          Type of boundary condition: 'v' or 'w'
%   c       [1 x 1]         Velocity of propagation

% OUTPUT
%   x       [2 x n(n-1)/2]  Coordinates of the mesh points x(:,i)=[x;t]
%   u       [2 x n(n-1)/2]  Solution in the mesh points u(:,i)=[v;w]


n = size(x_in, 2); % number of initial condition points

xb = x_in(1,1); % x coordinate of the left boundary (is it right to put 
% the bc in the first point of the intial condition points?)

x = zeros(2,n*(n-1)/2);
u = zeros(2,n*(n-1)/2);

k = 0;

if type == 'v'
    
    for i = 2:n
        x(1,i-1+k) = xb;
        x(2,i-1+k) = x_in(2,i) - ( xb - x_in(1,i) )/c;
    
        u(1,i-1+k) = ub; %v
        u(2,i-1+k) = ub - u_in(1,i) + u_in(2,i); %w
    
        for j = 1:n-i
            [x(:,i+k), u(:,i+k)] = innerProblem(x(:,i-1+k), x_in(:,i+j), ...
                                                u(:,i-1+k), u_in(:,i+j), c);
            k = k + 1;
        
        end
    
    end
    
elseif type == 'w'
    
    for i = 2:n
        x(1,i-1+k) = xb;
        x(2,i-1+k) = x_in(2,i) - ( xb - x_in(1,i) )/c;
    
        u(2,i-1+k) = ub; %w
        u(1,i-1+k) = ub - u_in(2,i) + u_in(1,i); %v
    
        for j = 1:n-i
            [x(:,i+k), u(:,i+k)] = innerProblem(x(:,i-1+k), x_in(:,i+j), ...
                                                u(:,i-1+k), u_in(:,i+j), c);
            k = k + 1;
        
        end
    
    end
else
    error('Wrong type of boundary condition')
end

end