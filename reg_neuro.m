
% Assignment 1 - Regularization

%% (a) Construct a 3D sphere of appropriate size

[X,Y,Z] = sphere; % Unit sphere presented with 20 tiles

%surf(X,Y,Z)
%axis equal

%% (b) Distribute 20 sources uniformly within the sphere 

N = 20; % Number of sources
eps = 0.1; % Distance to the boundary

% Calculate points
rvals = 2*rand(N,1)-1;
el = asin(rvals);
az = 2*pi*rand(N,1);
ra = (1-eps)*(rand(N,1).^(1/3));

[px,py,pz] = sph2cart(az,el,ra);

% Visualize

% In 3D
figure('Name',sprintf('3D-plot with eps = %0.2f',eps))
mesh(X,Y,Z,'FaceAlpha',0.5)
hold on
plot3(px,py,pz,'o')
hold off
axis('equal')
legend('Sphere','Sources')

% In 2D
figure('Name',sprintf('2D-plot with eps = %0.2f',eps))
% x-y
subplot(1,3,1)
plotcircle(1,0,0)
hold on
plot(px,py,'x')
hold off
axis('equal')
title('x-y-plane')
legend('Sphere','Sources','Location','South')
% x-z
subplot(1,3,2)
plotcircle(1,0,0)
hold on
plot(px,pz,'x')
hold off 
axis('equal')
title('x-z-plane')
legend('Sphere','Sources','Location','South')
% y-z
subplot(1,3,3)
plotcircle(1,0,0)
hold on
plot(py,pz,'x')
hold off
axis('equal')
title('y-z-plane')
legend('Sphere','Sources','Location','South')

%% (c) Simulate EEG recordings

tN = 101;

% Initiate variables 
Y1 = zeros(1,tN);
Y2 = zeros(1,tN);
Y3 = zeros(1,tN);
Y4 = zeros(1,tN);
Y5 = zeros(1,tN);

% Set Initial Conditions = 0

% Generate noise
w1 = randn(1,tN);
w2 = randn(1,tN);
w3 = randn(1,tN);
w4 = randn(1,tN);
w5 = randn(1,tN);

% Simulate data
for t = 4:tN
    
    Y1(t) = 0.95*sqrt(2)*Y1(t-1) - 0.9025*Y1(t-2) + w1(t);
    Y2(t) = 0.5*Y1(t-2) + w2(t);
    Y3(t) = -0.4*Y1(t-3) + w3(t);
    Y4(t) = -0.5*Y1(t-2) + 0.25*sqrt(2)*Y4(t-1) + w4(t);
    Y5(t) = 0.25*sqrt(2)*Y4(t-1) + 0.25*sqrt(2)*Y5(t-1) + w5(t);
    
end

%% (d) Distribute 5 electrodes at the surface of the sphere

M = 5; % Number of electrodes

% Calculate points
TH = 2*pi*rand(1,M);
PH = asin(-1+2*rand(1,M));
[pex,pey,pez] = sph2cart(TH,PH,1);

% Assign electrodes 1-5 to signals Y1-Y5
E1 = {pex(1), pey(1), pez(1), Y1};
E2 = {pex(2), pey(2), pez(2), Y2};
E3 = {pex(3), pey(3), pez(3), Y3};
E4 = {pex(4), pey(4), pez(4), Y4};
E5 = {pex(5), pey(5), pez(5), Y5};

% Visualize

% In 3D

figure('Name',sprintf('3D-plot with eps = %0.2f',eps))
mesh(X,Y,Z,'FaceAlpha',0.5)
hold on
plot3(px,py,pz,'o')
plot3(pex,pey,pez,'s','Color','r')
hold off
axis('equal')
legend('Sphere','Sources','Electrodes')

% In 2D
figure('Name',sprintf('2D-plot with eps = %0.2f',eps))
% x-y
subplot(1,3,1)
plotcircle(1,0,0)
hold on
plot(px,py,'x')
plot(pex,pey,'s','Color','r')
hold off
axis('equal')
title('x-y-plane')
legend('Sphere','Sources','Electrodes','Location','South')
% x-z
subplot(1,3,2)
plotcircle(1,0,0)
hold on
plot(px,pz,'x')
plot(pex,pez,'s','Color','r')
hold off 
axis('equal')
title('x-z-plane')
legend('Sphere','Sources','Electrodes','Location','South')
% y-z
subplot(1,3,3)
plotcircle(1,0,0)
hold on
plot(py,pz,'x')
plot(pey,pez,'s','Color','r')
hold off
axis('equal')
title('y-z-plane')
legend('Sphere','Sources','Electrodes','Location','South')

%% (e) Construct the lead field (gain matrix)

c = 2; % Constant c
G = zeros(M,N); % Initiate empty matrix

% Fill matrix
for m = 1:M
    
    for n = 1:N
        
        R = sqrt((pex(m)-px(n))^2+(pey(m)-py(n))^2+(pez(m)-pz(n))^2);
        G(m,n) = c / (R^2);
        
    end
    
end

%% (f) Solve the regularization problem for determining the intensities

x_star = zeros(tN,N); % Initiate empty x_star

% Fill x_star
for t = 1:tN
    
    V = [ Y1(t); Y2(t); Y3(t); Y4(t); Y5(t) ];
    x_star(t,:) = G.' * inv(G*G.') * V;
    
end

%% (g) Visualize the results

% Visualize solution in time of each source
% Distribute results to four different figures
for j = 1:4

    figure('Name',sprintf('Results: Part %d',j))
    for k = ((5*j)-4):(5*j)

        z = mod(k+4,5)+1;
        
        subplot(5,1,z)
        plot(0:tN-1,x_star(:,k))
        axis([0 tN-1 -inf inf])
        title(sprintf('Source: x=%0.4f, y=%0.4f, z=%0.4f', px(k),py(k),pz(k)))

    end

end

% Visualize signals of the electrodes
figure('Name','Signals of Electrodes')
subplot(5,1,1)
plot(0:tN-1,E1{4})
[phii,thetaa] = cart2pol(E1{1},E1{2},E1{3});
set(gca,'xticklabel',[])
title(sprintf('Position: \\theta =%0.4f, \\phi =%0.4f', phii, thetaa))
subplot(5,1,2)
plot(0:tN-1,E2{4})
[phii,thetaa] = cart2pol(E2{1},E2{2},E2{3});
set(gca,'xticklabel',[])
title(sprintf('Position: \\theta=%0.4f, \\phi=%0.4f', phii, thetaa))
subplot(5,1,3)
plot(0:tN-1,E3{4})
[phii,thetaa] = cart2pol(E3{1},E3{2},E3{3});
set(gca,'xticklabel',[])
title(sprintf('Position: \\theta=%0.4f, \\phi=%0.4f', phii, thetaa))
subplot(5,1,4)
plot(0:tN-1,E4{4})
[phii,thetaa] = cart2pol(E4{1},E4{2},E4{3});
set(gca,'xticklabel',[])
title(sprintf('Position: \\theta=%0.4f, \\phi=%0.4f', phii, thetaa))
subplot(5,1,5)
plot(0:tN-1,E5{4})
[phii,thetaa] = cart2pol(E5{1},E5{2},E5{3});
title(sprintf('Position: \\theta=%0.4f, \\phi=%0.4f', phii, thetaa))


