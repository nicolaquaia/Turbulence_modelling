close all
clc

% note: to run this you need the following files:
% chan180_means.txt, chan395_means.txt, chan590_means.txt
% chan180_reystress.txt, chan395_reystress.txt, chan590_reystress.txt
% folder named 'pictures' in the same directory

fontSize = 16;
filePath = 'pictures';

% variables

k = 0.41 ; % costant
rho = 1000 ; % water density
mu = 1e-3 ; % dynamic viscosity
nu = mu/rho ; % cinematic viscosity
H = 1 ; % channel's height
A_plus = 26;
B = 5.2;

n = 100;
Re_list = [180,395,590];


%% comparison velocity of DNS, Prandlt, damped function

path_list = {'chan180_means.txt', 'chan395_means.txt', 'chan590_means.txt'};
Re_list = [180,395,590];

% comparison for each Re_tau
for i = 1:length(Re_list)
    [y_ref, u_plus_ref, y_plus_ref] = read_profile_Umean(path_list{i});
    [u_plus_pra, y_pra, y_plus_pra] = compute_u_norm(Re_list(i), 1000, H, k);
    [u_plus_damp, y_damp, y_plus_damp] = compute_u_damp(Re_list(i), 1000, H, k, A_plus);

    figure('Units', 'centimeters', 'Position', [1, 1, 25, 10])
    hold on
    plot(y_plus_ref, u_plus_ref,y_plus_pra, u_plus_pra,y_plus_damp, u_plus_damp)
    legend('DNS','Prandtl', 'Damped', 'Location', 'southeast')
    set(gca, 'XScale', 'log', 'FontSize', fontSize);
    xlim([1, 1000])
    ylim([0, 25])
    xline(5, '--', 'Color', 'k', 'HandleVisibility', 'off');  % Dashed black line at y = 5
    xline(30, '--', 'Color', 'k', 'HandleVisibility', 'off'); % Dashed black line at y = 30
    text(12, 0.9 * max(ylim), mean([min(u_plus_ref), max(u_plus_ref)]), 'Buffer Layer', 'FontSize', 12, 'HorizontalAlignment', 'center');
    text(2, 0.9 * max(ylim), mean([min(u_plus_ref), max(u_plus_ref)]), 'Viscous sub-layer', 'FontSize', 12, 'HorizontalAlignment', 'center');
    text(80, 0.9 * max(ylim), mean([min(u_plus_ref), max(u_plus_ref)]), 'Log-law region', 'FontSize', 12, 'HorizontalAlignment', 'center');
    title(sprintf('Mean velocity profile for Re_{\\tau} = %d', Re_list(i)));
    ylabel('U+', 'FontSize', fontSize)
    xlabel('y+', 'FontSize', fontSize)
    grid on;
    hold off

    fileName = sprintf('comparison_%d.pdf', Re_list(i));

    exportgraphics(gcf, fullfile(filePath, fileName), 'ContentType', 'vector', 'BackgroundColor', 'none');
end


%% damping function analyisis

n=1000;
y_plus = linspace(0,200,n);
lt = k*y_plus;
lt_damp = zeros(n,1);
for i=1:n
    lt_damp(i) = k * y_plus(i) * (1 - exp(-y_plus(i)/A_plus));
end

figure('Units', 'centimeters', 'Position', [1, 1, 25, 10])
plot(y_plus, lt, y_plus, lt_damp);
ylabel('l_t', 'FontSize', fontSize)
xlabel('y+', 'FontSize', fontSize)
legend('simple', 'with damping', 'Location', 'southeast')
set(gca, 'FontSize', fontSize);
grid on
exportgraphics(gcf, fullfile(filePath, 'damping_analysis.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');


%% Reynold stresses


path_list = {'chan180_reystress.txt', 'chan395_reystress.txt', 'chan590_reystress.txt'};
[y_ref, y_plus_ref, R_uv_ref] = read_profile_reynold('chan590_reystress.txt');
[y, y_plus, uv] = compute_reynold(590, 1000, H, k);
color_list = {'r', 'b', 'g'};


figure('Units', 'centimeters', 'Position', [1, 1, 25, 10])
subplot(1, 2, 1);  % 1 row, 2 columns, first plot
hold on
for i = 1:3
    [y_ref, y_plus_ref, R_uv_ref] = read_profile_reynold(path_list{i});
    [y, y_plus, uv] = compute_reynold(Re_list(i), 1000, H, k);

    color = color_list{i};
    plot(y_plus_ref, -R_uv_ref,'--', 'Color', color)  % Dashed line for reference data
    plot(y_plus, -uv, 'Color', color)                  % Solid line for other data
end
hold off
xlim([0,80])
ylabel('R_{uv}', 'FontSize', fontSize)
xlabel('y+', 'FontSize', fontSize)
title('(a) Simple function')
set(gca,'FontSize', fontSize);
grid on

subplot(1, 2, 2);  % 1 row, 2 columns, first plot
hold on
for i = 1:3
    [y_ref, y_plus_ref, R_uv_ref] = read_profile_reynold(path_list{i});
    [y, y_plus, uv] = compute_reynold_damp(Re_list(i), 1000, H, k, A_plus);

    color = color_list{i};
    plot(y_plus_ref, -R_uv_ref,'--', 'Color', color)  % Dashed line for reference data
    plot(y_plus, -uv, 'Color', color)                % Solid line for other data
end
hold off

xlim([0,80])
xlabel('y+', 'FontSize', fontSize)
title('(b) With dampening')
grid on;
set(gca, 'FontSize', fontSize);
legend( '', 'Re_{\tau} = 180', '', 'Re_{\tau} = 395', '', 'Re_{\tau} = 590', 'Location', 'southeast')

exportgraphics(gcf, fullfile(filePath, 'Re_stress.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');



%% miscellaneous

n_list = [10, 100, 1000, 10000];

figure('Units', 'centimeters', 'Position', [1, 1, 25, 10])
hold on
for i = 1:length(n_list)
    [u_plus, y, y_plus] = compute_u_norm(180, n_list(i), H, k);
    plot(y_plus, u_plus, 'DisplayName', sprintf('n = %d', n_list(i)))
end
hold off
set(gca, 'XScale', 'log', 'FontSize', fontSize);
ylabel('U+', 'FontSize', fontSize)
xlabel('y+', 'FontSize', fontSize)
grid on;
legend('Location', 'northwest'); 
exportgraphics(gcf, fullfile(filePath, 'various_n.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');


figure('Units', 'centimeters', 'Position', [1, 1, 25, 10])
hold on
for i = 1:length(Re_list)
    [u_plus, y, y_plus] = compute_u_norm(Re_list(i), 1000, H, k);
    plot(y_plus, u_plus, 'DisplayName', sprintf('Re_{\\tau} = %d', Re_list(i)))
end
hold off
set(gca, 'XScale', 'log', 'FontSize', fontSize);
ylabel('U+', 'FontSize', fontSize)
xlabel('y+', 'FontSize', fontSize)
grid on;
legend('Location', 'northwest'); 
exportgraphics(gcf, fullfile(filePath, 'various_Re.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');


%% 

[y, u_ref, y_plus] = read_profile_Umean('chan590_means.txt');
n = length(y_plus);

u_near = y_plus;

u_log = zeros(n,1);
for i=1:n
    u_log(i) = 1/k* log(y_plus(i)) + B;
end


figure('Units', 'centimeters', 'Position', [1, 1, 25, 10])
hold on
plot(y_plus, u_ref)
plot(y_plus, u_near)
plot(y_plus, u_log)
hold off
ylim([0,30])
%xlim([1,50])
set(gca, 'XScale', 'log', 'FontSize', fontSize);
ylabel('U+', 'FontSize', fontSize)
xlabel('y+', 'FontSize', fontSize)
legend('DNS', 'u+=y+', 'log law', 'Location', 'northwest')
exportgraphics(gcf, fullfile(filePath, 'log_law.pdf'), 'ContentType', 'vector', 'BackgroundColor', 'none');


%% comparison with laminar

Re_list = [180,395,590];
comp_coeff = 5;

for i=1:3
    
    Re_tau = Re_list(i);
    u_tau = Re_tau*nu/H;
    [y_lam, u_lam] = compute_u_lam(n, H, nu, u_tau);
    [y_turb, u_turb] = compute_u_turb(n, H, k,nu,u_tau);

    figure('Units', 'centimeters', 'Position', [1, 1, 25, 10])
    hold on
    plot(u_lam, y_lam, 'LineWidth',1.5)
    plot(comp_coeff*u_turb,y_turb,'LineWidth',1.5)
    set(gca, 'FontSize', fontSize);
    xlabel('U [m/s]','FontSize', fontSize)
    ylabel('y [m]','FontSize', fontSize)
    legend('Laminar','Turbulent','Location','southeast')
    grid minor
    title(sprintf('Re_{\\tau} = %d', Re_list(i)));


    fileName = sprintf('laminar_%d.pdf', Re_list(i));
    exportgraphics(gcf, fullfile(filePath, fileName), 'ContentType', 'vector', 'BackgroundColor', 'none');
end


%% functions
function [A_mat] = A_forward(n)
    uni = ones(1,n); % vector of ones

    A = ones(n);
    A = diag(-1*ones(1,n)) + diag(uni(1:end-1), -1) ;
    A(1,2) = 0;
    A(end,end-1) = 0 ;
    A(1,1) = 1 ;
    A(end,end) = 1 ;
    A(end,end-1) = -1 ;

    A_mat = A;
    %disp(A_mat)
end

function [u_plus, y, y_plus] = compute_u_norm(Re_tau, n, H, k)

    y = linspace(0,H,n);
    y_plus = y * Re_tau/H;
    dy_plus = diff(y_plus);
    dy_plus = dy_plus(1);

    A = A_forward(n);

    b = zeros(n,1);
    for i = 1:n
        b(i) = (-1 + sqrt(1 + 4 * k^2 * y_plus(i)^2 * (1 - y_plus(i)/Re_tau))) / (2 * k^2*y_plus(i)^2);
    end
    b(1) = 0;
    b(end) = 0;
    b = b*dy_plus;

    u_plus = -A\b;
end


function [u_plus, y, y_plus] = compute_u_damp(Re_tau, n, H, k, A_plus)

    y = linspace(0,H,n);
    y_plus = y * Re_tau/H;
    dy_plus = diff(y_plus);
    dy_plus = dy_plus(1);

    A = A_forward(n);

    b = zeros(n,1);
    for i = 1:n
        damp = 1 - exp(-y_plus(i)/A_plus);
        b(i) = (-1 + sqrt(1 + 4 * k^2 * y_plus(i)^2 *damp^2 * (1 - y_plus(i)/Re_tau))) / (2 * k^2*y_plus(i)^2 * damp^2);
    end
    b(1) = 0;
    b(end) = 0;
    b = b*dy_plus;

    u_plus = -A\b;
end


function [y, y_plus, uv] = compute_reynold(Re_tau, n, H, k)

    y = linspace(0,H,n);
    y_plus = y * Re_tau/H;

    dudy = zeros(n,1);
    for i = 1:n
        dudy(i) = (-1 + sqrt(1 + 4 * k^2 * y_plus(i)^2 * (1 - y_plus(i)/Re_tau))) / (2 * k^2*y_plus(i)^2);
    end

    uv = dudy' + (y_plus/Re_tau) - 1;
end


function [y, y_plus, uv] = compute_reynold_damp(Re_tau, n, H, k, A_plus)

    y = linspace(0,H,n);
    y_plus = y * Re_tau/H;

    dudy = zeros(n,1);
    for i = 1:n
        damp = 1 - exp(-y_plus(i)/A_plus);
        dudy(i) = (-1 + sqrt(1 + 4 * k^2 * y_plus(i)^2 *damp^2* (1 - y_plus(i)/Re_tau))) / (2 * k^2*y_plus(i)^2*damp^2);
    end

    uv = dudy' + (y_plus/Re_tau) - 1;
end

function [y,u] = compute_u_lam(n, H, nu, u_tau)
    y = linspace(0,H,n);
    u = zeros(n,1);

    for i=1:n
        u(i) = -((u_tau^2)/(nu*H))*(0.5*y(i)^2) + y(i)*((u_tau^2) / nu);
    end
end


function [y,u] = compute_u_turb(n, H, k,nu,u_tau)
    y = linspace(0,H,n);
    dy = diff(y);
    dy = dy(1);

    A = A_forward(n);

    b = zeros(n,1);
    for i = 1:n
        lt = k * y(i) ;
        b(i) = (-nu + sqrt(nu^2 + 4 * lt^2 * u_tau^2 * (1-(y(i)/(H))))) / (2*lt^2) ;
    end
    b(1) = 0;
    b(end) = 0;
    b = b*dy;

    u = -A\b;
end


function [y, Umean, y_plus] = read_profile_Umean(path) 
    file = fopen(path, 'r');
    % Skip the first lines
    for i = 1:25
        fgetl(file);
    end
    
    % Read the data
    data = textscan(file, '%f %f %f %f %f %f %f');
    
    % Close the file
    fclose(file);
    
    % Access columns
    y = data{1};
    Umean = data{3};
    y_plus = data{2};
end


function [y, y_plus, R_uv] = read_profile_reynold(path) 
    file = fopen(path, 'r');
    % Skip the first lines
    for i = 1:25
        fgetl(file);
    end
    
    % Read the data
    data = textscan(file, '%f %f %f %f %f %f %f %f %f');
    
    % Close the file
    fclose(file);
    
    % Access columns
    y = data{1};
    y_plus = data{2};
    R_uv = data{6};
end