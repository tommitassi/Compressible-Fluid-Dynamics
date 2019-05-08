%% Method of Characteristics - Acoustics 1D

%% 5. Creation of the upper mesh

function [ x, u ] = upperMesh( x_l, u_l, x_r, u_r, c )

% upperMesh

% INPUT
%   x_l     [2 x nn]        Coordinates of the left boundary mesh points x(:,i)=[x;t]
%   u_l     [2 x nn]        Solution in the left boundary mesh points u(:,i)=[v;w]
%   x_r     [2 x nn]        Coordinates of the right boundary mesh points x(:,i)=[x;t]
%   u_r     [2 x nn]        Solution in the right boundary mesh points u(:,i)=[v;w]
%   c       [1 x 1]         Velocity of propagation

% OUTPUT
%   x       [2 x N]         Coordinates of the upper mesh points x(:,i)=[x;t]
%   u       [2 x N]         Solution in the upper mesh points u(:,i)=[v;w]

nn = size(x_l, 2);
n = (-1+sqrt(1+8*nn))/2;
N = n*n;

x = zeros(2, N);
u = zeros(2, N);

l = 1;
p = 1;

for i = 1 : n
    
    q = 1;
    
    for j = 1 : n
        
        [x(:,l), u(:,l)] = innerProblem( x_l(:,p), x_r(:,q),...
                                         u_l(:,p), u_r(:,q), c);
                                         
        l = l + 1;
        
        q = q + n-j+1;
                                           
    end
    
    p = p + n-i+1;
        
end


end