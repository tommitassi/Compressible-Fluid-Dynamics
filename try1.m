clear

c = 340;

x_inp = 0:50:2000;
n = length(x_inp);

x_in = [x_inp; zeros(1,n)];

u_in = zeros(2,n);



[ x_un, u_un ] = unboundedDomain( x_in, u_in, c);
[ x_l, u_l ] = leftBoundary_Moving( x_in, u_in, 1, 'v', 500, 100, c );
[ x_r, u_r ] = rightBoundary( x_in, u_in, 1, 'w', c );
%[ x_up, u_up ] = upperMesh( x_l, u_l, x_r, u_r, c );

close all
figure
plot(x_un(1,:),x_un(2,:),'o')
hold on
plot(x_l(1,:),x_l(2,:),'o')
plot(x_r(1,:),x_r(2,:),'o')
%plot(x_up(1,:),x_up(2,:),'o')
% quiver(x_un(1,:),x_un(2,:),u_un(1,:),u_un(2,:))

