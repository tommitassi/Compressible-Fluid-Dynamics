%% Method of Characteristics - Acoustics 1D

%% 1. Inner problem

function [ I, u ] = innerProblem( Pl, Pr, ul, ur, c )

% innerProblem finds the intesection point of the characteristic lines
% starting from two points and evaluate the solution in that point

% INPUT
%   Pl      [2 x 1]         Coordinates of the left point (P+) [x; t]
%   Pr      [2 x 1]         Coordinates of the right point (P-) [x; t]
%   ul      [2 x 1]         Solution in the left point [v; w]
%   ur      [2 x 1]         Solution in the right point [v; w]
%   c       [1 x 1]         Velocity of propagation

% OUTPUT
%   I       [2 x 1]         Coordinates of the intersection [x; t]
%   u       [2 x 1]         Solution in the intersection [v; w]

I = zeros(2,1);

I(1) = ( Pr(1) + Pl(1) )/2 - c*( Pl(2) - Pr(2) )/2;
I(2) = ( Pr(1) - Pl(1) )/( 2 * c ) + ( Pr(2) + Pl(2) )/2;

u = zeros(2,1);

u(1) = ( ur(1) + ul(1) )/2 + ( ul(2) - ur(2) )/2;
u(2) = ( ur(2) + ul(2) )/2 + ( ul(1) - ur(1) )/2;

end