%% Method of Characteristics - Acoustics 1D

%% 4. Creation of the mesh relative to the left moving boundary:

function [ x, u ] = leftmoving_Boundary( x_in, u_in, ub, type, x_f, v0, c )

% leftmoving_Boundary permits to complete the mesh using the left
% boundary condition and to compute the solution in those points.
% The boundary moves with a constant deceleration until the braking
% distance x_f, from which it remains still.

% INPUT
%   x_in    [2 x n]         Coordinates of the points where the initial
%                           conditions are imposed x_in(:,i) = [x;t]
%   u_in    [2 x n]         Initial conditions u_in(:,i) = [v;w]
%   ub      [1 x 1]         Boundary condition
%   type    string          Type of boundary condition: 'v' or 'w'
%   x_f     [1 x 1]         x-braking distance from the initial point [m]
%   v0      [1 x 1]         Initial velocity of the boundary [m/s]
%   c       [1 x 1]         Velocity of propagation [m/s]

% OUTPUT
%   x       [2 x n(n-1)/2]  Coordinates of the mesh points x(:,i) = [x;t]
%   u       [2 x n(n-1)/2]  Solution in the mesh points u(:,i) = [v;w]

% Parabola: x = A*x^2 + v0*x:
A = - v0^2 / (4 * x_f);

% Braking time:
t_f = (-v0 + sqrt(v0^2 + 4*A*x_f))/(2*A); % [s]

% Space coordinate from which the still boundary is used:
x_lim = x_f + c*t_f; % [m]

% Number of initial condition points:
n = size(x_in, 2); 

% Initialization of two vectors:
x = zeros(2,n*(n-1)/2);
u = zeros(2,n*(n-1)/2);

k = 0;

if type == 'v'
    
    for i = 2:n
        if x_in(1,i) <= x_lim
            x(2,i-1+k) = ( -c-v0+sqrt((v0+c)^2+4*A*(c*x_in(2,i)+x_in(1,i))) )/(2*A);
            x(1,i-1+k) = x_in(1,i) - c*(x(2,i-1+k)-x_in(2,i));
            
        else
            x(1,i-1+k) = x_f;
            x(2,i-1+k) = x_in(2,i) - ( x_f - x_in(1,i) )/c;
    
        end
            
        u(1,i-1+k) = ub; % v
        u(2,i-1+k) = ub - u_in(1,i) + u_in(2,i); % w
            
        for j = 1:n-i
            [x(:,i+k), u(:,i+k)] = innerProblem(x(:,i-1+k), x_in(:,i+j), ...
                                                u(:,i-1+k), u_in(:,i+j), c);
            k = k + 1;                            
            
        end
            
    end
    
elseif type == 'w'
    
    for i = 2:n
        if x_in(1,i) <= x_lim
            x(2,i-1+k) = ( -c-v0+sqrt((v0+c)^2+4*A*(c*x_in(2,i)+x_in(1,i))) )/(2*A);
            x(1,i-1+k) = x_in(1,i) - c*(x(2,i-1+k)-x_in(2,i));
            
        else
            x(1,i-1+k) = x_f;
            x(2,i-1+k) = x_in(2,i) - ( x_f - x_in(1,i) )/c;
            
        end
    
        u(2,i-1+k) = ub; % w
        u(1,i-1+k) = ub - u_in(2,i) + u_in(1,i); % v
            
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