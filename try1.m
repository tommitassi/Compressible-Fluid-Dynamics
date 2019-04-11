clear all

c = 2;

x_inp = [1 2 3 4];

x_in = [x_inp; 1 1 1 1 ];

u_in = [1 2 3 4; 2 2 2 2];

n = size(u_in,2);

[ x_un, u_un ] = unboundedDomain( x_in, u_in, c);
[ xl, ul ] = leftBoundary( x_in, u_in, 1, 'w', c );


close all
figure
plot(x_un(1,:),x_un(2,:),'o')
hold on
plot(xl(1,:),xl(2,:),'o')
% quiver(x_un(1,:),x_un(2,:),u_un(1,:),u_un(2,:))

