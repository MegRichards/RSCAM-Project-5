clear all
clc
%% create polygon
poly1 = polyshape([0 0 2 2],[2 0 0 2]);
poly2 = polyshape([1 2 2 1],[2 2 1 1]);
polygon = subtract(poly1,poly2);
%polygon = polyshape([0 0 1 1],[1 0 0 1]); % rectangle

% center and rescale the polygon
[x_c,y_c] = centroid(polygon);
polygon = translate(polygon,-[x_c,y_c]);
corners = polygon.Vertices;
scl = max([ max(corners(:,1))-min(corners(:,1)) , max(corners(:,2))-min(corners(:,2)) ]);
polygon = scale(polygon,1/scl);
%% Implement correct pole placement

% establish constants outside the loop
f = @(z) real(z).^2; % function for boundary data
z_star = 0;
sigma = 4;
eps = 1e-3; % Desired accuracy of solution
noise = 0.3; % level of randomness for poles
distmin = 1e-15; % min distance from pole to corner
distmax = 1e-1; % max distance from pole to new boundary point

% get input values
corners = polygon.Vertices;
number_corners = length(polygon.Vertices);

% Increase # of poles at corners iteratively
Error = [];
% create n new poles at each iteration
% for increasing n with sqrt(n) approx evenly spaced
% n_vals = [25, 81, 169, 289, 441, 625, 841, 1089];
n_vals = linspace(1,8,10).^2;
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
                poles = [poles new_pole(polygon, ncorner, beta, noise*rand, false)];
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
                poles = [poles new_pole(polygon, ncorner, beta, noise*rand, false)];
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
