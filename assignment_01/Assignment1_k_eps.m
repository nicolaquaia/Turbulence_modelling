close all
clc

fontSize = 16;
filePath = 'pictures';

n = 1000;
Re_tau = 180;

H=1;

C_mu = 0.09;
C1 = 1.45;
C2 = 1.9;
sigma_k = 1;
sigma_e = 1.3;
A_plus = 26;

rho = 1000 ; % water density
mu = 1e-3 ; % dynamic viscosity
nu = mu/rho ; % cinematic viscosity
%nu = 3.5000e-4 ; % cinematic viscosity

i_max = 100;
err_min = 1e-10;
alpha = 0.1;
beta = 0.4;
check = 0;
ind = 0;
rel_start = 4;

[y, y_plus, u, u_plus] = k_eps(n, H, Re_tau, nu, C_mu,C1,C2,sigma_k,sigma_e,A_plus, ...
    i_max,err_min,alpha,beta,check,ind,rel_start);

[y_d, y_plus_d, u_d, u_plus_d] = k_eps_damp(n, H, Re_tau, nu, C_mu,C1,C2,sigma_k,sigma_e,A_plus, ...
    i_max,err_min,alpha,beta,check,ind,rel_start);


figure('Units', 'centimeters', 'Position', [1, 1, 25, 10])
hold on
subplot(1, 2, 1); 
plot(u, y,u_d, y_d);
hold on
xlabel('U [m/s]', 'FontSize', fontSize);
ylabel('y [m]', 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
grid on;


subplot(1, 2, 2); 
hold on
semilogx(y_plus, u_plus, y_plus_d, u_plus_d);
hold off
xlabel('y+','FontSize', fontSize);
ylabel('U+','FontSize', fontSize);
set(gca, 'XScale', 'log', 'FontSize', fontSize);
legend('k-eps', 'damped k-eps','Location', 'northwest')
grid on;

exportgraphics(gcf, fullfile(filePath, 'k_eps.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');

% zooom

figure('Units', 'centimeters', 'Position', [1, 1, 25, 10])
hold on
subplot(1, 2, 1); 
plot(u, y,u_d, y_d);
hold off
ylim([0.9, 1])
xlabel('U [m/s]', 'FontSize', fontSize);
ylabel('y [m]', 'FontSize', fontSize);
set(gca, 'FontSize', fontSize);
grid on;

subplot(1, 2, 2); 
hold on
semilogx(y_plus, u_plus, y_plus_d, u_plus_d);
hold off
xlim([160, 180])
xlabel('y+','FontSize', fontSize);
ylabel('U+','FontSize', fontSize);
set(gca, 'XScale', 'log', 'FontSize', fontSize);
legend('k-eps', 'damped k-eps','Location', 'northwest')
grid on;

exportgraphics(gcf, fullfile(filePath, 'k_eps_zoom.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');





%% functions


function [nu_t] = get_nu_t(n, nu_t, k, dis, C_mu)
    for i = 1:n
        nu_t(i) = C_mu * k(i)^2 / dis(i);
    end
end


function [u_out] = get_u(n, u, nu_t, dy, nu, Re_tau, H, alpha, beta, ind, rel_start)
    a1 = zeros(n,1);
    a2 = zeros(n,1);
    a3 = zeros(n,1);
    b = zeros(n,1);
    A = zeros(n);
    
    for i = 2:n-1
        a1(i) = 2*nu + nu_t(i-1) + nu_t(i);
        a2(i) = -(4*nu + nu_t(i+1) + 2*nu_t(i) + nu_t(i-1));
        a3(i) = 2*nu + nu_t(i+1) + nu_t(i);
        b(i) = -2*dy^2*nu^2 * Re_tau^2 / H^3;
    end
    
    %[b, a2] = relaxation_before(n, b, a2, u, alpha,ind, rel_start);
    
    %disp(a1)
    for j=2:n-1
        A(j,j) = a2(j);
        A(j,j-1) = a1(j);
        A(j,j+1) = a3(j);
    end

    % boundaries
    A(1,1) = 1;
    A(end,end) = 1;
    A(end,end-1) = -1;

    b(1) = 0;
    b(end) = 0;

    u_new = A\b;
    u_out = relaxation_after(u, u_new, beta);
end


function [prod] = get_prod(n, prod, u, nu_t, dy)
    for i = 2:n-1
        prod(i) = nu_t(i) * ( (u(i+1) - u(i-1))/(2*dy) )^2;
    end

    % boundaries ???

end


function [k_out] = get_k(n, k, nu_t, dis, dy, nu, sigma_k, alpha, beta, ind, rel_start)
    a1 = zeros(n,1);
    a2 = zeros(n,1);
    a3 = zeros(n,1);
    b = zeros(n,1);
    A = zeros(n);
    
    for i = 2:n-1
        a1(i) = sigma_k*2*nu + nu_t(i-1) + nu_t(i);
        a2(i) = -(sigma_k*4*nu + nu_t(i+1) + 2*nu_t(i) + nu_t(i-1));
        a3(i) = sigma_k*2*nu + nu_t(i+1) + nu_t(i);
        b(i) = 2*dy^2*sigma_k*(dis(i) - prod(i));
    end

    %[b, a2] = relaxation_before(n, b, a2, k, alpha,ind, rel_start);
    
    %disp(a1)
    for j=2:n-1
        A(j,j) = a2(j);
        A(j,j-1) = a1(j);
        A(j,j+1) = a3(j);
    end

    % boundaries
    A(1,1) = 1;
    A(end,end) = 1;
    A(end,end-1) = -1;

    b(1) = 0;
    b(end) = 0;

    k_new = A\b;
    k_out = relaxation_after(k, k_new, beta);
end

function [dis_out] = get_dis(n, dis, nu_t, k, dy, nu, sigma_e, C1, C2, alpha, beta, ind, rel_start)
    a1 = zeros(n,1);
    a2 = zeros(n,1);
    a3 = zeros(n,1);
    b = zeros(n,1);
    A = zeros(n);
    
    for i = 2:n-1
        a1(i) = sigma_e*2*nu + nu_t(i-1) + nu_t(i);
        a2(i) = -(sigma_e*4*nu + nu_t(i+1) + 2*nu_t(i) + nu_t(i-1));
        a3(i) = sigma_e*2*nu + nu_t(i+1) + nu_t(i);
        b(i) = 2*dy^2*sigma_e* (C2*dis(i)^2/k(i) - C1*dis(i)/k(i)*prod(i));
    end

    %[b, a2] = relaxation_before(n, b, a2, dis, alpha,ind, rel_start);
    
    for j=2:n-1
        A(j,j) = a2(j);
        A(j,j-1) = a1(j);
        A(j,j+1) = a3(j);
    end

    % boundaries
    A(1,1) = 1;
    A(end,end) = 1;
    A(end,end-1) = -1;

    b(1) = 2*nu*k(1)/(dy^2);            % why?
    %b(end) = 2*nu*k(end-1)/(dy^2);      % why?
    b(end) = 0;

    dis_new = A\b;
    dis_out = relaxation_after(dis, dis_new, beta);
end


function [u] = get_u_poiseuille(n, dy, nu, Re_tau, H)
    a1 = zeros(n,1);
    a2 = zeros(n,1);
    a3 = zeros(n,1);
    b = zeros(n,1);
    A = zeros(n);
    
    for i = 2:n-1
        a1(i) = 1;
        a2(i) = -2;
        a3(i) = 1;
        b(i) = - nu * Re_tau^2 * dy^2 / (H^3);
    end
    
    for j=2:n-1
        A(j,j) = a2(j);
        A(j,j-1) = a1(j);
        A(j,j+1) = a3(j);
    end

    % boundaries
    A(1,1) = 1;
    A(end,end) = 1;
    A(end,end-1) = -1;

    b(1) = 0;
    b(end) = 0;

    u = A\b;
end


function [res] = relaxation_after(x, x_new, beta)
    res = beta * x_new + (1-beta)*x;
end


function [b, a2] = relaxation_before(n, b, a2, x, alpha,ind,rel_start)
    if ind > rel_start
        b = b + a2.*x.*(1-alpha)/alpha;
        a2 = a2./alpha;
    end

    % for i =2:n-1
    %     b(i) = b(i) + a2(i)*x(i)*(1-alpha)/alpha;
    %     a2(i) = a2(i)/alpha;
    % end
end



function [y, y_plus, u, u_plus] = k_eps(n, H, Re_tau, nu, C_mu,C1,C2,sigma_k,sigma_e,A_plus, ...
    i_max,err_min,alpha,beta,check,ind,rel_start)

    u_tau = Re_tau*nu/H;
    y = linspace(0,H,n);
    y_plus = y*u_tau/nu;
    dy = diff(y);
    dy = dy(1);
    
    % initialize arrays
    nu_t = zeros(n,1);
    u = zeros(n,1);         % initial guess can be improved
    prod = zeros(n,1);
    k = ones(n,1);
    dis = ones(n,1);
    f_mu = zeros(n,1);

    while check==0
        
        ind = ind+1;
    
        % compute nu_t
        nu_t_new = get_nu_t(n, nu_t, k, dis, C_mu);
        
        % compute u
        u_new = get_u(n, u, nu_t, dy, nu, Re_tau, H, alpha, beta, ind, rel_start);
        % compute prod
        prod_new = get_prod(n, prod, u, nu_t, dy);
        % compute k
        k_new = get_k(n, k, nu_t, dis, dy, nu, sigma_k, alpha, beta, ind, rel_start);
        % compute dis
        dis_new = get_dis(n, dis, nu_t, k, dy, nu, sigma_e, C1, C2, alpha, beta, ind, rel_start);
    
    
        % compute error
        err = sum((u-u_new).^2);
    
        % update
        nu_t = nu_t_new;
        u = u_new;
        prod = prod_new;
        k = k_new;
        dis = dis_new;
    
        % check
        if err < err_min
            check = 1;
        end
        if ind > i_max-1
            check = 1;
        end
       
        %disp('end loop')
        fprintf('i: %5d, Error: %15.10f\n', ind, err);
    end

    u_plus = u/u_tau;

end



function [y, y_plus, u, u_plus] = k_eps_damp(n, H, Re_tau, nu, C_mu,C1,C2,sigma_k,sigma_e,A_plus, ...
    i_max,err_min,alpha,beta,check,ind,rel_start)

    u_tau = Re_tau*nu/H;
    y = linspace(0,H,n);
    y_plus = y*u_tau/nu;
    dy = diff(y);
    dy = dy(1);
    
    % initialize arrays
    nu_t = zeros(n,1);
    u = zeros(n,1);         % initial guess can be improved
    prod = zeros(n,1);
    k = ones(n,1);
    dis = ones(n,1);
    f_mu = zeros(n,1);

    while check==0
        
        ind = ind+1;
    
        % compute nu_t
        nu_t_new = get_nu_t(n, nu_t, k, dis, C_mu);
    
        % damping function
        for i=1:n
            f_mu(i) = (1 - exp(-y_plus(i)/A_plus))^2;
            nu_t_new(i) = f_mu(i)*nu_t_new(i);
        end
    
        % compute u
        u_new = get_u(n, u, nu_t, dy, nu, Re_tau, H, alpha, beta, ind, rel_start);
        % compute prod
        prod_new = get_prod(n, prod, u, nu_t, dy);
        % compute k
        k_new = get_k(n, k, nu_t, dis, dy, nu, sigma_k, alpha, beta, ind, rel_start);
        % compute dis
        dis_new = get_dis(n, dis, nu_t, k, dy, nu, sigma_e, C1, C2, alpha, beta, ind, rel_start);
    
    
        % compute error
        err = sum((u-u_new).^2);
    
        % update
        nu_t = nu_t_new;
        u = u_new;
        prod = prod_new;
        k = k_new;
        dis = dis_new;
    
        % check
        if err < err_min
            check = 1;
        end
        if ind > i_max-1
            check = 1;
        end
       
        %disp('end loop')
        fprintf('i: %5d, Error: %15.10f\n', ind, err);
    end

    u_plus = u/u_tau;

end