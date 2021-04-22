function [weights, err, A] = solve_LSQR(type1, type2, varargin)

if type1 == "Laplace"
    if type2 == "interior"
        % eg solve_LSQR(square, Laplace, interior, [p1..],z_star, N2, f = @(x) 2*x^3+7*x^2+x, [bp1..]);
        poles = varargin{1}; % poles are outside domain (row vector, important!)
        d = varargin{2}; % distances of poles from their respective corner (row vector, important!)
        z_star = varargin{3}; % point at centre of domain 
        N2 = varargin{4}; % max degree of polynomail terms
        f = varargin{5}; % function for Dirichlet boundary data
        bps = varargin{6}; % these bps are boundary points (ie zk) (column vector, important!)
        b = f(bps);
        A =@(p) [real(d./(p-poles)) imag(d./(p-poles)) real((p-z_star).^(0:N2)) imag((p-z_star).^(1:N2))];
        A_ = A(bps);
        weights = (A_'*A_)\(A_'*b);
        err = max(abs(b-A_*weights))/max(abs(b));
    
    elseif type2 == "exterior"
        poles = varargin{1}; % poles are inside domain (row vector, important!)
        d = varargin{2}; % distances of poles from their respective corner (row vector, important!)
        z_star = varargin{3}; % point at centre of domain 
        N2 = varargin{4}; % max degree of polynomail terms
        f = varargin{5}; % function for Dirichlet boundary data
        bps = varargin{6}; % these bps are boundary points (ie zk) (column vector, important!)
        b = -f(bps);
        scl = min(abs(bps-z_star)); % scaling parameter for numerical stability
        A =@(p) [real(d./(p-poles)) imag(d./(p-poles)) real((scl./(p-z_star)).^(1:N2)) imag((scl./(p-z_star)).^(1:N2))];
        A_ = A(bps);
        weights = (A_'*A_)\(A_'*b);
        err = max(abs(b-A_*weights))/max(abs(b));
        
    elseif type2 == "transmission"
        % arguments for interior
        poles_in = varargin{1}; % poles are outside domain (row vector, important!)
        d_in = varargin{2}; % distances of poles for int problem from their respective corner (row vector, important!)
        % arguments for exterior
        poles_ex = varargin{3}; % poles are inside domain (row vector, important!)
        d_ex = varargin{4}; % distances of poles for ext problem from their respective corner (row vector, important!)
        H = varargin{5}; % harmonic function that u tends to at infinity
        grad_H = varargin{6}; % gradient of harmonic ft
        % joint arguments
        N2 = varargin{7}; % max degree of polynomail terms
        bps = varargin{8}; % these bps are boundary points (ie zk) (column vector, important!)
        norms = varargin{9}; % normals associated with bps (column vector, important!)
        z_star = varargin{10}; % point at centre of domain 
        kappa = varargin{11}; % interior kappa value in problem, set exterior kappa to 1
        % get value for A by combining int and ex probs and solving on
        % Get b
        dH_dn = real(norms).*real(grad_H(bps)) + imag(norms).*imag(grad_H(bps)); 
        b = [H(bps); dH_dn];
        % boundary condition 1
        A_in =@(p) [real(d_in./(p-poles_in)) imag(d_in./(p-poles_in)) real((p-z_star).^(0:N2)) imag((p-z_star).^(1:N2))];
        A_in_ = A_in(bps);
        scl=1;
        A_ex =@(p) [real(d_ex./(p-poles_ex)) imag(d_ex./(p-poles_ex)) real((scl./(p-z_star)).^(1:N2)) imag((scl./(p-z_star)).^(1:N2))];
        A_ex_ = A_ex(bps);
        % boundary condition 2
        r_i_real = [real(-d_in./(bps-poles_in).^2) imag(-d_in./(bps-poles_in).^2) real((0:N2).*((bps-z_star).^(-1:(N2-1)))) imag((1:N2).*((bps-z_star).^(0:(N2-1))))];
        r_i_imag = [imag(-d_in./(bps-poles_in).^2) -real(-d_in./(bps-poles_in).^2) imag((0:N2).*((bps-z_star).^(-1:(N2-1)))) -real((1:N2).*((bps-z_star).^(0:(N2-1))))];
        A_i_diff = real(norms).*r_i_real - imag(norms).*r_i_imag;
        
        r_e_real = [real(-d_ex./(bps-poles_ex).^2) imag(-d_ex./(bps-poles_ex).^2) real(-(1:N2).*((scl./(bps-z_star)).^(2:(N2+1)))) imag(-(1:N2).*((scl./(bps-z_star)).^(2:(N2+1))))];
        r_e_imag = [imag(-d_ex./(bps-poles_ex).^2) -real(-d_ex./(bps-poles_ex).^2) imag(-(1:N2).*((scl./(bps-z_star)).^(2:(N2+1)))) -real(-(1:N2).*((scl./(bps-z_star)).^(2:(N2+1))))];
        A_e_diff = real(norms).*r_e_real - imag(norms).*r_e_imag;
        
        % Combine
        A_ = [A_in_ -A_ex_;kappa*A_i_diff -A_e_diff];
        
        % outputs
        A= {A_in,A_ex}; 
        weights = (A_'*A_)\(A_'*b);
        %weights = lsqr(A_,b,1e-8,1000);
        err = max(abs(b-A_*weights))/max(abs(b));
        weights = {weights(1:size(A_in_,2)), weights((size(A_in_,2)+1):end)};
    end
    
elseif type1 == "Helmholtz"
    if type2 == "interior"
    elseif type2 == "exterior"
        poles = varargin{1}; % poles are outside domain (row vector, important!)
        d = varargin{2}; % distances of poles from their respective corner (row vector, important!)
        z_star = varargin{3}; % point at centre of domain 
        N2 = varargin{4}; % max degree of polynomail terms
        f = varargin{5}; % function for Dirichlet boundary data
        bps = varargin{6}; % these bps are boundary points (ie zk) (column vector, important!)
        k = varargin{7};
        b = -f(bps);
        A =@(p)  [d.*besselh(-1,k*abs(p-poles)).*((p-poles)./abs(p-poles)).^(-1) d.*besselh(1,k*abs(p-poles)).*((p-poles)./abs(p-poles))...
            H_j(N2,k*abs(p-z_star)).*((p-z_star)./abs(p-z_star)).^(-N2:N2)];
        A_ = A(bps);
        weights = (A_'*A_)\(A_'*b);
        err = max(abs(b-A_*weights))/max(abs(b));
    end
end
