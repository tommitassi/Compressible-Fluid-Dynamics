%% Method of Characteristics - Acoustics 1D

%% 1. Definition of the physical problem and results:

clear
close all

% Physical data:
c = 340; % Speed of sound [m/s]
v0 = 100; % Initial velocity of the train [m/s]
x_f = 200; % Braking distance [m]

x_inp = 0:50:2000; 
n = length(x_inp); 

% Space-time coordinates of the initial data line:
x_in = [x_inp; zeros(1,n)]; % [x;t]

% Initial condition along the initial data line:
u_in = zeros(2,n); % [v;w]

% Procedure:
[ x_un, u_un ] = unbounded_Domain(x_in, u_in, c);
[ x_l, u_l ] = leftmoving_Boundary(x_in, u_in, 10, 'v', x_f, v0, c);
[ x_r, u_r ] = right_Boundary(x_in, u_in, 1, 'w', c);
[ x_up, u_up ] = upper_Mesh( x_l, u_l, x_r, u_r, c );

% Graphic representation of the results:
figure
plot(x_un(1,:),x_un(2,:),'o')
hold on
plot(x_l(1,:),x_l(2,:),'o')
plot(x_r(1,:),x_r(2,:),'o')
plot(x_up(1,:),x_up(2,:),'o')

% COMMENTI E NOTE SUL LAVORO MANCANTE:

% 1. Costruzione della soluzione per tempi successivi in corrispondenza del
%    left boundary.

% 2. Risolvere il 'problema dell'ascoltatore': data una determinata
%    posizione dello spazio valutare l'evoluzione della perturbazione nel
%    tempo generando un grafico perturbazione-tempo a spazio fissato.

% 3. Rappresentazione grafica della soluzione ottenuta in un grafico
%    tridimensionale (provare ad utilizzare la funzione stem3).