clear all
clc
%% create polygon
poly1 = polyshape([0 0 2 2],[2 0 0 2]);
poly2 = polyshape([1 2 2 1],[2 2 1 1]);
polygon = subtract(poly1,poly2);
polygon = polyshape([0 0 1 1],[1 0 0 1]); % rectangle
n = 6; polygon = polyshape(real(exp(1i*[0:n-1]/n*2*pi))',imag(exp(1i*[0:n-1]/n*2*pi))'); % 

% center and rescale the polygon 
[x_c,y_c] = centroid(polygon);
polygon = translate(polygon,-[x_c,y_c]);
corners = polygon.Vertices;
scl = max([ max(corners(:,1))-min(corners(:,1)) , max(corners(:,2))-min(corners(:,2)) ]);
polygon = scale(polygon,1/scl);
%% Interior Dirichlet Laplace
% establish constants outside the loop
f = @(z) real(z).^2.*imag(z); % function for boundary data
z_star = 0;
sigma = 4;
noise = 0.5; % level of randomness for poles
distmin = 1e-13; % min distance from pole to corner
distmax = 1e-1; % max distance from pole to new boundary point
% get input values
corners = polygon.Vertices;
number_corners = length(polygon.Vertices);

% Increase # of poles at corners iteratively
Error = [];
% create n new poles at each iteration
% for increasing n with sqrt(n) approx evenly spaced
% n_vals = [25, 81, 169, 289, 441, 625, 841, 1089];
n_vals = linspace(2,10,5).^2;
% iterate
for n = n_vals
    poles = [];
    d = []; % distances of poles from their respective corner
    
    normal = [];
    bps = []; % boundary points
    [bps, normal] = discretization_boundary(polygon, 50); %IS THIS ALWAYS 30?
    for ncorner=1:number_corners
        % get interior angle
        angle = angle_check(polygon,ncorner);
        %CHECK STRICT INEQ HERE, might be easier to use phi (new_pole_ang)
        if angle < 0
            % add 3n poles
            for j=1:(3*n)
                beta = exp(-sigma*(sqrt(3*n)-sqrt(j)));
                if beta>distmin
                d = [d beta];
                poles = [poles new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), false)];
                end
                if beta<distmax
                    delta = beta;
                    [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, delta/3);
                    [b3, b4, n3, n4] = new_boundary_normal(polygon, ncorner, (2*delta)/3);
                    [b5, b6, n5, n6] = new_boundary_normal(polygon, ncorner, delta);
                    bps = [bps ; b1 ; b2; b3; b4; b5; b6];
                end
            end
        else
            % at each iteration, add n poles contained in beta_val
            for j=1:n
                beta = exp(-sigma*(sqrt(n)-sqrt(j))); % beta calc from rule 3.2 in ref 2
                if beta>distmin
                d = [d beta];
                poles = [poles new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), false)];
                end
                if beta<distmax
                    delta = beta;
                    [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, delta/3);
                    [b3, b4, n3, n4] = new_boundary_normal(polygon, ncorner, (2*delta)/3);
                    [b5, b6, n5, n6] = new_boundary_normal(polygon, ncorner, delta);
                    bps = [bps ; b1 ; b2; b3; b4; b5; b6];
                end
            end
        end
    end
    disp(n)
    N2 = n;
    [~,err,~] = solve_LSQR('Laplace', 'interior', poles, d, z_star, N2, f, bps);
    Error = [Error err];
end
figure, plot(Error), set(gca, 'YScale', 'log') % plot error

% Plot the final solution
[weights,~,A] = solve_LSQR('Laplace', 'interior', poles, d, z_star, N2, f, bps);
plot_solution(polygon, 'Laplace', 'interior', A, weights)
plot(real(poles),imag(poles),'.')
axis equal;

%% Exterior Dirichlet Laplace
% establish constants outside the loop
f = @(z) real(z); % function for boundary data
z_star = 0;
sigma = 4;
noise = 0; % level of randomness for poles
distmin = 1e-15; % min distance from pole to corner
distmax = 1e-1; % max distance from pole to new boundary point

% get input values
corners = polygon.Vertices;
number_corners = length(polygon.Vertices);

% Increase # of poles at corners iteratively
Error = [];
% create n new poles at each iteration
% for increasing n with sqrt(n) approx evenly spaced
% n_vals = [25, 81];
n_vals = linspace(2,10,5).^2;
% iterate
for n = n_vals
    poles = [];
    d = []; % distances of poles from their respective corner
    
    normal = [];
    bps = []; % boundary points
    [bps, normal] = discretization_boundary(polygon, 50); %IS THIS ALWAYS 30?
    for ncorner=1:number_corners
        % get interior angle
        angle = angle_check(polygon,ncorner);
        %CHECK STRICT INEQ HERE, might be easier to use phi (new_pole_ang)
        if angle < 0
            % add 3n poles
            for j=1:(3*n)
                beta = exp(-sigma*(sqrt(3*n)-sqrt(j)));
                if (beta>distmin && beta<0.5)
                d = [d beta];
                poles = [poles new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), true)];
                end
                if beta<distmax
                    delta = beta;
                    [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, delta/3);
                    [b3, b4, n3, n4] = new_boundary_normal(polygon, ncorner, (2*delta)/3);
                    [b5, b6, n5, n6] = new_boundary_normal(polygon, ncorner, delta);
                    bps = [bps ; b1 ; b2; b3; b4; b5; b6];
                end
            end
        else
            % at each iteration, add n poles contained in beta_val
            for j=1:n
                beta = exp(-sigma*(sqrt(n)-sqrt(j))); % beta calc from rule 3.2 in ref 2
                if (beta>distmin && beta<0.5)
                d = [d beta];
                poles = [poles new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), true)];
                end
                if beta<distmax
                    delta = beta;
                    [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, delta/3);
                    [b3, b4, n3, n4] = new_boundary_normal(polygon, ncorner, (2*delta)/3);
                    [b5, b6, n5, n6] = new_boundary_normal(polygon, ncorner, delta);
                    bps = [bps ; b1 ; b2; b3; b4; b5; b6];
                end
            end
        end
    end
    disp(n)
    N2 = round(n)+1;
    if length(poles)
    [~,err,~] = solve_LSQR('Laplace', 'exterior', poles, d, z_star, N2, f, bps);
    end
    Error = [Error err];
end    

figure, plot(Error), set(gca, 'YScale', 'log') % plot error    
    
% Plot the final solution
[weights,~,A] = solve_LSQR('Laplace', 'exterior', poles, d, z_star, N2, f, bps);
plot_solution(polygon, 'Laplace', 'exterior', A, weights, f), hold on
plot(real(poles),imag(poles),'.')
axis equal;

%% Exterior Dirichlet Helmholtz
% establish constants outside the loop
k = 30; % wavenumber
%source = 1/2 + 1*1i; f = @(z) besselh(0,k*abs(z-source)); % point source
dir = [1 -1]; dir  = dir/norm(dir); f = @(z) exp(1i*k*(dir(1)*real(z)+dir(2)*imag(z))); % plane wave

% establish constants outside the loop
z_star = 0;
sigma = 4;
eps = 1e-3; % Desired accuracy of solution
noise = 0.; % level of randomness for poles
distmin = 1e-15; % min distance from pole to corner
distmax = 2e-1; % max distance from pole to new boundary point

% get input values
corners = polygon.Vertices;
number_corners = length(polygon.Vertices);

% Increase # of poles at corners iteratively
Error = [];
% create n new poles at each iteration
% for increasing n with sqrt(n) approx evenly spaced
%n_vals = [2, 8, 16, 28, 44, 62];
n_vals = linspace(2,10,5).^2;
% iterate
for n = n_vals
    poles = [];
    d = []; % distances of poles from their respective corner
    %z_star = [z_star 0.5*(2*(rand-0.5)+1i*2*(rand-0.5))];
    normal = [];
    bps = []; % boundary points
    [bps, normal] = discretization_boundary(polygon, 50); %IS THIS ALWAYS 30?
    for ncorner=1:number_corners
        % get interior angle
        angle = angle_check(polygon,ncorner);
        %CHECK STRICT INEQ HERE, might be easier to use phi (new_pole_ang)
        if angle < 0
            % add 3n poles
            for j=1:(3*n)
                beta = exp(-sigma*(sqrt(3*n)-sqrt(j)));
                if (beta>distmin && beta<0.5)
                d = [d beta];
                poles = [poles new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), true)];
                end
                if beta<distmax
                    delta = beta;
                    [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, delta/3);
                    [b3, b4, n3, n4] = new_boundary_normal(polygon, ncorner, (2*delta)/3);
                    [b5, b6, n5, n6] = new_boundary_normal(polygon, ncorner, delta);
                    bps = [bps ; b1 ; b2; b3; b4; b5; b6];
                end
            end
        else
            % at each iteration, add n poles contained in beta_val
            for j=1:n
                beta = exp(-sigma*(sqrt(n)-sqrt(j))); % beta calc from rule 3.2 in ref 2
                if (beta>distmin && beta<0.5)
                d = [d beta];
                poles = [poles new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), true)];
                end
                if beta<distmax
                    delta = beta;
                    [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, delta/3);
                    [b3, b4, n3, n4] = new_boundary_normal(polygon, ncorner, (2*delta)/3);
                    [b5, b6, n5, n6] = new_boundary_normal(polygon, ncorner, delta);
                    bps = [bps ; b1 ; b2; b3; b4; b5; b6];
                end
            end
        end
    end
    disp(n)
    N2 = round(n/2)+1;
    if length(poles)
    [~,err,~] = solve_LSQR('Helmholtz', 'exterior', poles, d, z_star, N2, f, bps, k);
    end
    Error = [Error err];
end    

figure, plot(Error), set(gca, 'YScale', 'log') % plot error    
    
% Plot the final solution
[weights,~,A] = solve_LSQR('Helmholtz', 'exterior', poles, d, z_star, N2, f, bps, k);
plot_solution(polygon, 'Helmholtz', 'exterior', A, weights, f), hold on
plot(real(poles),-imag(poles),'.');
axis equal;
%% Interior Dirichlet Laplace
% establish constants outside the loop
f = @(z) real(z).^2.*imag(z); % function for boundary data
z_star = 0;
sigma = 4;
noise = 0; % level of randomness for poles
distmin = 1e-13; % min distance from pole to corner
distmax = 1e-1; % max distance from pole to new boundary point
% get input values
corners = polygon.Vertices;
number_corners = length(polygon.Vertices);

% Increase # of poles at corners iteratively
Error = [];
% create n new poles at each iteration
% for increasing n with sqrt(n) approx evenly spaced
% n_vals = [25, 81, 169, 289, 441, 625, 841, 1089];
n_vals = linspace(2,5,5).^2;
% iterate
for n = n_vals
    poles = [];
    d = []; % distances of poles from their respective corner
    
    normal = [];
    bps = []; % boundary points
    [bps, normal] = discretization_boundary(polygon, 100); %IS THIS ALWAYS 30?
    for ncorner=1:number_corners
        % get interior angle
        angle = angle_check(polygon,ncorner);
        %CHECK STRICT INEQ HERE, might be easier to use phi (new_pole_ang)
        if angle < 0
            % add 3n poles
            for j=1:(3*n)
                beta = exp(-sigma*(sqrt(3*n)-sqrt(j)));
                if (beta>distmin && beta<0.5)
                d = [d beta];
                poles = [poles new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), false)];
                end
                if beta<distmax
                    delta = beta;
                    [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, delta/3);
                    [b3, b4, n3, n4] = new_boundary_normal(polygon, ncorner, (2*delta)/3);
                    [b5, b6, n5, n6] = new_boundary_normal(polygon, ncorner, delta);
                    bps = [bps ; b1 ; b2; b3; b4; b5; b6];
                end
            end
        else
            % at each iteration, add n poles contained in beta_val
            for j=1:n
                beta = exp(-sigma*(sqrt(n)-sqrt(j))); % beta calc from rule 3.2 in ref 2
                if (beta>distmin && beta<0.5)
                d = [d beta];
                poles = [poles new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), false)];
                end
                if beta<distmax
                    delta = beta;
                    [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, delta/3);
                    [b3, b4, n3, n4] = new_boundary_normal(polygon, ncorner, (2*delta)/3);
                    [b5, b6, n5, n6] = new_boundary_normal(polygon, ncorner, delta);
                    bps = [bps ; b1 ; b2; b3; b4; b5; b6];
                end
            end
        end
    end
    disp(n)
    N2 = n;
    [~,err,~] = solve_LSQR('Laplace', 'interior', poles, d, z_star, N2, f, bps);
    Error = [Error err];
end
figure, plot(Error), set(gca, 'YScale', 'log') % plot error

% Plot the final solution
[weights,~,A] = solve_LSQR('Laplace', 'interior', poles, d, z_star, N2, f, bps);
plot_solution(polygon, 'Laplace', 'interior', A, weights)
plot(real(poles),imag(poles),'.')
axis equal;
%% Transmission
% Get arguments constant across interior and exterior 
z_star = 0;
sigma = 4;
noise = 0; % level of randomness for poles
distmin = 1e-13; % min distance from pole to corner
distmax = 1e-2; % max distance from pole to new boundary point
%H = @(z) real(z.^2); % function for boundary data
%grad_H = @(z) real(2*z)-1i*imag(2*z);
H = @(z) real(z);
grad_H = @(z) 1;
kappa = 0.1;  % not sure what value kappa should be?
% get input values
corners = polygon.Vertices;
number_corners = length(polygon.Vertices);

% Increase # of poles at corners iteratively
Error = [];
n_vals = linspace(2,8,5).^2;
% iterate
for n = n_vals
    poles_in = [];
    d_in = []; % distances of poles from their respective corner
    poles_ex = [];
    d_ex = [];
    
    normal = [];
    bps = []; % boundary points
    [bps, normal] = discretization_boundary(polygon, 100); %IS THIS ALWAYS 30?
    for ncorner=1:number_corners
        % get interior angle
        angle = angle_check(polygon,ncorner);
        %CHECK STRICT INEQ HERE, might be easier to use phi (new_pole_ang)
        if angle < 0
            % add 3n poles
            for j=1:(3*n)
                beta = exp(-sigma*(sqrt(3*n)-sqrt(j)));
                if (beta>distmin && beta<.5)
                d_in = [d_in beta];
                poles_in = [poles_in new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), false)];
                d_ex = [d_ex beta];
                poles_ex = [poles_ex new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), true)];
                end
                if beta<distmax
                    delta = beta;
                    [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, delta/3);
                    [b3, b4, n3, n4] = new_boundary_normal(polygon, ncorner, (2*delta)/3);
                    [b5, b6, n5, n6] = new_boundary_normal(polygon, ncorner, delta);
                    bps = [bps ; b1 ; b2; b3; b4; b5; b6];
                    normal = [normal; n1 ; n2; n3; n4; n5; n6];
                end
            end
        else
            % at each iteration, add n poles contained in beta_val
            for j=1:n
                beta = exp(-sigma*(sqrt(n)-sqrt(j))); % beta calc from rule 3.2 in ref 2
                if (beta>distmin && beta<.5)
                d_in = [d_in beta];
                poles_in = [poles_in new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), false)];
                d_ex = [d_ex beta];
                poles_ex = [poles_ex new_pole(polygon, ncorner, beta, noise*2*(rand-0.5), true)];
                end
                if beta<distmax
                    delta = beta;
                    [b1, b2, n1, n2] = new_boundary_normal(polygon, ncorner, delta/3);
                    [b3, b4, n3, n4] = new_boundary_normal(polygon, ncorner, (2*delta)/3);
                    [b5, b6, n5, n6] = new_boundary_normal(polygon, ncorner, delta);
                    bps = [bps ; b1 ; b2; b3; b4; b5; b6];
                    normal = [normal; n1 ; n2; n3; n4; n5; n6];
                end
            end
        end
    end
    N2 = round(n)+1;
    disp(N2);
    if length(poles_in)
    [~,err,~] = solve_LSQR_trans('Laplace', 'transmission', poles_in, d_in, poles_ex,d_ex,H,grad_H, N2, bps, normal, z_star, kappa);
    Error = [Error err];
    end
end    

figure, plot(Error), set(gca, 'YScale', 'log') % plot error  

% Plot the final solution
[weights,~,A] = solve_LSQR_trans('Laplace', 'transmission', poles_in, d_in, poles_ex,d_ex,H,grad_H, N2, bps, normal, z_star, kappa);
plot_solution(polygon, 'Laplace', 'transmission', A, weights, H)
plot(real(poles_in),imag(poles_in),'.', 'color', 'blue'), hold on
plot(real(poles_ex),imag(poles_ex),'.', 'color', 'red')
axis equal
