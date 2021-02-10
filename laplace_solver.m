% Laplace equation solver on polygonal domains
% Arguments:
%   w = List of boundary vertices, counterclockwise
%   center = Point centrally situated within domain used for interior
%   polynomial-based approximation
%   h = Boundary value problem (Dirichlet), specified per-edge as functions
%   of the points on the boundary
%   tol = Desired accuracy of solution
% Return Value:
%   u = Solution to boundary value problem, valid on the domain only
% Optional Named Parameter Pairs:
%   'tests' = false, Turn on visual inspection of boundary and poles before
%   computations.
%   'neumann' = false, Turn on support for neumann boundary conditions.
%   Note, this does not work with solely neumann, at least one must be
%   dirichlet. Specify a side is neumann by adding 1 to second column, 0 if
%   dirichlet, i.e. [{@(z) z} 1; repmat([{@(z) 0} 0]), length(w), 1]]
%   'discont' = false, Weights approximation and error by distance to 
%   nearest corner, allowing for convergence with corner discontinuities.
%   'curved' = false, Allows for curved boundaries, specified in curves.
%   Does not support neumann boundaries at same time currently.
%   'curves' = [], Parameterisation of curved boundaries. Example row would
%   contain [{@(t) exp(1i*(t+pi/4))} 1.5*pi 1], the function, maximal value
%   for t and the normal at the corner.
%   'curved_hull' = [], A polygonal representation of the boundary for use
%   in inclusion testing and plotting.
%   'plot3' = false, Plot using contours or 3D representation.

function u = laplace_solver(w, center, h, tol, varargin)
    p = inputParser; % Handle optional parameters
    addOptional(p, 'tests', false);
    addOptional(p, 'neumann', false); % Does not support curves yet
    addOptional(p, 'discont', false);
    addOptional(p, 'curved', false)
    addOptional(p, 'curves', []);
    addOptional(p, 'curved_hull', []);
    addOptional(p, 'plot3', false);
    addOptional(p, 'heatmap', false);
    addOptional(p, 'solfunc', false);
    addOptional(p, 'plots', true);
    addOptional(p, 'fixed', false);
    addOptional(p, 'sigma', 4);
    
    parse(p, varargin{:});
    neumann = p.Results.neumann;
    curved = p.Results.curved;
    if curved, curves = p.Results.curves; end
    if curved && neumann, disp("Error: Curved domains are not supported with neumann boundaries yet"); end
    
    % Get the maximal diameter of the domain
    w_r = sort(real(w)); w_r = w_r([1 end]);
    w_i = sort(imag(w)); w_i = w_i([1 end]);
    scale = max([diff(w_r),diff(w_i)]);

    if curved % Get hull-based scale
        b_r = sort(real(p.Results.curved_hull)); b_r = b_r([1 end]);
        b_i = sort(imag(p.Results.curved_hull)); b_i = b_i([1 end]);
        scale = max([diff(b_r),diff(b_i)]);
    end

    % Confirm validity of poles and boundary sampling points if desired
    if p.Results.tests
        clf
        poles = compute_poles_test(w, scale, 3, center, p);
        y = input('next? ');
        bdd = bdd_sample_test(w, poles, h, 3, tol, center, p); 
        y = input('next? ');
    end

    % Increase n until convergence to specified tolerance.
    max_steps = 40;
    prev_err = inf; bdd = []; flat_poles = [];
    errs = inf(1,max_steps);conds = zeros(1,max_steps); Nvec = zeros(1,max_steps);
    for n = 3:max_steps
        % Exp. clustered poles at corners
        poles = compute_poles(w, scale, n, p);
        flat_poles = cell2mat(poles).';
        [~, N1] = size(flat_poles);

        % Distance of poles from corners and boundary sample points
        [bdd, b, dist] = bdd_sample(w, poles, h, n, tol, p);
        [M,~] = size(bdd);
        [Ms, ~] = size(b);
        newman = dist./(bdd - flat_poles);
        
        % Polynomial approx in bulk
        N2 = ceil(n/2);
        runge = ((bdd - center)./scale).^(1:N2);

        % Construct least squares matrix
        N = 2*N1 + 2*N2;
        disp("n = "+n); disp("N = "+N); disp("N1 = "+N1); disp("N2 = "+N2);
        disp("M = "+M); disp("Ms = "+Ms);

        A = [real(runge) imag(runge) real(newman) imag(newman) ones(M,1)];

        % If necessary, weight by distance to vertex
        weights = ones(1,M);
        if p.Results.discont, weights = min(abs(bdd.' - w)); end

        % Handle Neumann Boundaries
        D = eye(M, M);
        b_neumann = zeros(M,1);
        if neumann
            D = zeros(Ms, M);
            zero_i = 1;
            for i = 1:Ms
                if b(i,2) == 1
                    D(i,:) = [zeros(1, zero_i-1) 1/tol -1/tol zeros(1, M-zero_i-1)];
                    b_neumann(zero_i+1) = 1;
                    zero_i = zero_i + 2;
                else
                    D(i,:) = [zeros(1, zero_i-1) 1 zeros(1, M-zero_i)];
                    zero_i = zero_i + 1;
                end
            end
        end
        weights = weights(b_neumann == 0);
        weight_m = spdiags(sqrt(weights.'), 0, Ms, Ms);
        
        % Solve using backslash, matricies are always ill-conditioned
        warn = warning('off','MATLAB:rankDeficientMatrix');
        
        % Truncated QR-Solver
        % Uses StrongRRQR, found at https://math.berkeley.edu/~mgu/MA273/Strong_RRQR.pdf
        % Implemented by https://www.mathworks.com/matlabcentral/fileexchange/69139-strong-rank-revealing-qr-decomposition
%         f = 1.03;
%         [Q, R, perm] = sRRQR(D*A, f, 'tol', 10^-12);
%         [~,perminv] = sort(perm);
%         c = (R\(Q.'*b(:,1)));
%         c = c(perminv);

        % Normal QR-Solver
%         [Q,R,perm] = qr(D*A);
%         c = perm*(R\(Q.'*b(:,1)));
        
        % Normal SVD Solver
%         [U,S,V] = svd(D*A);
%         Truncate it 
%         truncSVD = @(U,S,V,p) V(:,1:p)*diag(1./diag(S(1:p,1:p)))*U(:,1:p)';
%         c = truncSVD(U,S,V,N)*b(:,1);
%         c = V*(S\U.'*b(:,1));
%         
        % Backslash solver
        c = (D * A) \ (b(:,1));
        
        warning(warn.state,'MATLAB:rankDeficientMatrix');
        % Check error on boundary
        err_v = weight_m*(D*A*c - b(:,1));
        err = norm(err_v,inf);
        disp("err="+err);
%         disp("cond="+cond(D*A));
%         conds(n) = cond(D*A);
        % Exit if tolerance reached or if error begins to increase
        % as convergence is not possible due to numerical instability.
        errs(n) = err; Nvec(n) = n;
        if ~p.Results.fixed && err <= tol %|| err > prev_err*10
            break;
        end
        prev_err = err;
    end
    errs = errs(Nvec > 0);
    conds = conds(Nvec > 0);
    Nvec = Nvec(Nvec > 0);
    
    % Construct proposed solution
    u = @(z) [real(((z-center)./scale).^(1:N2)) ...
            imag(((z-center)./scale).^(1:N2))   ...
            real(dist./(z-flat_poles))        ...
            imag(dist./(z-flat_poles)) ones(size(z))]*c;
    
    % Plot solution, poles and boundary sampling points
    % Plotting Config
    LW = 'linewidth'; MS = 'markersize'; FS = 'fontsize';
    fs = 16; PO = 'position'; FW = 'fontweight'; NO = 'normal';
    lw=1;
    if curved, w_r = b_r; w_i = b_i; w = p.Results.curved_hull; end
    ax = [w_r(1:2); w_i(1:2)] + .2*scale*[-1 1 -1 1]';
    axwide = [w_r(1:2); w_i(1:2)] + 1.1*scale*[-1 1 -1 1]';

    % Create sampling grid
    sx = linspace(w_r(1),w_r(2),100); sy = linspace(w_i(1),w_i(2),100);
    [xx,yy] = meshgrid(sx,sy); zz = xx + 1i*yy;

    % Evaluate solution on grid, discard points not in domain
    ff = zz; ff(:) = u(zz(:)); ff(~indomain(zz,w,p)) = NaN;
    
    % Get bounds for color scale, plot solution, domain and poles/bdd
    clf, shg
    figsize = [0,0,10,8];
    fontsize = 16;
    h = figure;
    
    levels = linspace(min(min(real(ff))),max(max(real(ff))),20);
    if p.Results.plots % Begin Plots
    
    if p.Results.tests
        contour(sx,sy,ff,levels,LW, 1), colorbar, hold on;
        plot(w([1:end 1]), '-k', LW,1), plot(flat_poles, '.r', MS,12), ...
            plot(bdd, '.g', MS, 12);
        plot(real(center),imag(center),'.k',MS,12), hold off;
        y = input('next? '); 

    end
    clf, shg
    axes(PO,[.4 .25 .6 .6])
    if p.Results.plot3, plot3(xx,yy,ff), colorbar;
    else
        contour(sx,sy,ff,levels,LW,.5), colorbar, axis equal, hold on;
        plot(w([1:end 1]), '-k', LW,1), plot(flat_poles, '.r', MS,10);
        set(gca,FS,fs-1,LW,lw), plot(real(center),imag(center),'.k',MS,12), axis(ax)
    end
    title(['#bdd sample points = ' int2str(M) ...
        ', #poles = ' int2str(N1)],FS,fs,FW,NO, ...
        "Interpreter", "latex"), hold off;

    set(gca, 'Units', 'normalized',...
        'FontUnits', 'points',...
        'FontWeight', 'normal',...
        'FontSize', fontsize, ...
        'FontName', 'Times',...
        'TickLabelInterpreter', 'latex',...
        LW,lw);
    
    % Plot boundary error
    axes(PO,[.05 .6 .25 .2])
    errmin = .01*tol; 
    ws = 'error'; if p.Results.discont, ws = 'weighted error'; end
    axis([-pi pi .0001*errmin 1]), grid on

    semilogy(angle(bdd(b_neumann == 0) - center), abs(err_v),'.k',MS,6), hold off

    set(gca,'ytick',10.^(-16:4:0))
    set(gca,'xtick',pi*(-1:1),'xticklabel',{'-\pi','0','\pi'})
    set(gca,FS,fs-1,LW,lw), xlabel('angle on boundary wrt $z_*$',FS,fs, ...
        "Interpreter", "latex"), ylabel(ws, FS, fs)
    title([ws ' on boundary'],FS,fs,FW,NO,"Interpreter", "latex")

    set(gca, 'Units', 'normalized',...
        'FontUnits', 'points',...
        'FontWeight', 'normal',...
        'FontSize', fontsize, ...
        'FontName', 'Times',...
        'TickLabelInterpreter', 'latex');
    % Plot speed of convergence
    axes(PO,[.05 .2 .25 .2]); ws = 'log-error'; if p.Results.discont, ws = 'weighted log-error'; end
    semilogy(sqrt(Nvec),errs,'.-k',LW,0.7,MS,10), grid on, hold on
%     semilogy(sqrt(Nvec),1./conds,'.-r',LW,0.7,MS,10), grid on, hold on
    errmin = .01*tol; axis([0.9*min(sqrt(Nvec)) 1.1*max(sqrt(Nvec)) errmin 10])
    set(gca,FS,fs-1,LW,lw), title('convergence',FS,fs,FW,NO, ...
        "Interpreter", "latex")
    xlabel('$\sqrt{n}$',FS,fs, ...
        "Interpreter", "latex"), ylabel(ws,FS,fs)
    set(gca,'ytick',10.^(-16:4:0),LW,lw)
    
    set(h, 'Units','Inches');
    pos = get(h,'Position');
    set(h,...
        'Position', figsize,...
        'Units', 'Inches', ...
        'PaperPositionMode','Auto',...
        'PaperUnits','Inches')
    set(gca, 'Units', 'normalized',...
        'FontUnits', 'points',...
        'FontWeight', 'normal',...
        'FontSize', fontsize, ...
        'FontName', 'Times',...
        'TickLabelInterpreter', 'latex');
    print(h,'poster\figures\curved-L','-depsc')

    end % End Plotting
    % TODO: A posteriori error on finer boundary sampling
    u = {sx,sy,ff, err, errs, n};
end

% Utility to check if point is in domain
function in = indomain(z, w, p)
    % TODO: make this more accurate for curved boundaries
    if p.Results.curved, w = p.Results.curved_hull; end
    in = inpolygon(real(z),imag(z),real(w),imag(w));
end

% Generate exponentially clustered poles along vertex normals of domain
% TODO: Vectorize computations further
function poles = compute_poles(w, scale, n, p)
    [m, ~] = size(w);
    sigma = p.Results.sigma;
    curved = p.Results.curved;
    if curved, curves = p.Results.curves; end

    % Edge normals
    n_e = zeros(m, 1);

    for i = 1:m
        n_e(i) = (w(i) - w(mod(i + m - 2, m)+1))*-1i;
        n_e(i) = n_e(i)/norm(n_e(i));
    end

    % Vertex normals (external angle bisectors)
    n_v = zeros(m, 1);

    for i = 1:m
        n_v(i) = n_e(i) + n_e(mod(i, m)+1);
        n_v(i) = n_v(i)/norm(n_v(i));
    end

    % Pole generation
    poles = cell(m,1);
    re_entrant = real(n_v).*real(n_e*1i)+imag(n_v).*imag(n_e*1i);
    if curved
        re_entrant = -ones(m,1);
        n_v = cell2mat(curves(:,3));
    end
    for i = 1:m
        N = n;
        % Check for re-entrancy
        if re_entrant(i) < 0
            N = 3*n;
        end

        dist = exp(-sigma*(sqrt(N) - sqrt(1:N)))';
        pole = w(i) + n_v(i) * scale .* dist;
        % Check for poles in domain
        poles{i} = pole(~indomain(pole,w,p));
    end
end

% Generate boundary sample points and get pole scale (distance to corners)
% TODO: Vectorize computations further
function [bdd,b,dist] = bdd_sample(w, poles, h, n, tol, p)
    [m, ~] = size(w);
    curved = p.Results.curved;
    neumann = p.Results.neumann;

    if curved, curves = p.Results.curves; end

    if neumann, n = n*2; end

    bdd = zeros(3*n*m,1);
    dist = zeros(3*n*m,1);
    b = zeros(3*n*m,2);
    bdd_i=1; dist_i=1; b_i = 1;

    for i = 1:m
        [N, ~] = size(poles{i});
        prev = mod(m + i - 2, m)+1; next = mod(i, m)+1;
        neg = w(prev) - w(i); pos = w(next) - w(i);

        % Setup for neumann bcs
        n_norm = imag(neg) - real(neg)*1i; n_norm = n_norm / abs(n_norm);
        p_norm = -imag(pos) + real(pos)*1i; p_norm = p_norm / abs(p_norm);
        if neumann
            if h{prev,2} == 1, L_n = 0:1; else, L_n = 0; end
            if h{i,2} == 1, L_p = 0:1; else, L_p = 0; end
        else
            L_n = 0; L_p = 0;
        end

        if curved, neg = curves{prev,2}; pos = curves{i,2}; end
        % loop through poles, checking distance
        for j = 1:N
            dist(dist_i) = norm(w(i) - poles{i}(j));
            for k = 1:3
                t = k/3 * dist(dist_i);
                eps = tol;
                if(t < norm(neg))
                    for l = L_n
                        bdd(bdd_i) = t * neg/norm(neg) + w(i) ...
                            + l*eps*n_norm;
                        if curved, bdd(bdd_i) = ...
                            curves{prev, 1}(curves{prev, 2} - t); end
                        if l == 0
                            b(b_i, 1) = h{prev}(bdd(bdd_i));
                            if neumann, b(b_i, 2) = h{prev,2}; end
                            b_i = b_i + 1;
                        end
                        bdd_i = bdd_i + 1;
                    end
                end
                if(t < norm(pos))
                    for l = L_p
                        bdd(bdd_i, 1) = t * pos/norm(pos) + w(i) ...
                           + l*eps*p_norm;
                        if curved, bdd(bdd_i, 1) = curves{i, 1}(t); end
                        if l == 0
                            b(b_i,1) = h{i}(bdd(bdd_i));
                            if neumann, b(b_i,2) = h{i,2}; end
                            b_i = b_i + 1;
                        end
                        bdd_i = bdd_i + 1;
                    end
                end
            end
            dist_i = dist_i+1;
        end
    end
    bdd = bdd(1:(bdd_i-1),:);
    dist = dist(1:(dist_i-1),:).';
    b = b(1:(b_i-1),:);
end

% Tests for the above functions
function out = compute_poles_test(w, scale, n, center, p)
    figsize = [0,0,2,2];
    fontsize = 9;
    h = figure;
    
    LW = 'linewidth'; MS = 'markersize';
    plot(w([1:end 1]), '-k', LW,1);
    out = compute_poles(w, scale, n, p);
    flat_poles = cell2mat(out);
    hold on, plot(flat_poles, '.r', MS,7), plot(real(center),imag(center),'.k',MS,7), hold off;
    ylim([-0.5 2.5])
    xlim([-0.5 2.5])
    set(h, 'Units','Inches');
        pos = get(h,'Position');
        set(h,...
            'Position', figsize,...
            'Units', 'Inches', ...
            'PaperPositionMode','Auto',...
            'PaperUnits','Inches')
        set(gca, 'Units', 'normalized',...
            'FontUnits', 'points',...
            'FontWeight', 'normal',...
            'FontSize', fontsize, ...
            'FontName', 'Times',...
            'TickLabelInterpreter', 'latex');
        print(h,'report\figures\domain','-depsc')
end

function out = bdd_sample_test(w, poles, h, n, tol, center, p)
    [out,b,dist] = bdd_sample(w, poles, h, n, tol, p);
    LW = 'linewidth'; MS = 'markersize';
    hold on, plot(out(:,1), '.g', MS,7); 
    plot(real(center),imag(center),'.k',MS,7), hold off;
end
