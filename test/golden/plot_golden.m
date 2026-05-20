% Plot rabe.nc output from the golden record run.
% Usage: plot_golden        (reads rabe.nc in current directory)
%        plot_golden myfile.nc

pkg load netcdf;

if nargin >= 1
    path = argv(){1};
else
    path = 'rabe.nc';
end

s             = ncread(path, 's_tor');
Lambda_A     = ncread(path, 'Lambda_A');
Lambda_B     = ncread(path, 'Lambda_B');
Lambda_S = ncread(path, 'Lambda_S');

info = ncinfo(path);
varnames = {info.Variables.Name};
has_sc = any(strcmp(varnames, 'lambda_LC_bB'));

if has_sc
    lambda_LC = ncread(path, 'lambda_LC_bB');
    remainder = ncread(path, 'remainder');
    figure('Position', [100 100 900 700]);
else
    figure('Position', [100 100 900 350]);
end

if has_sc
    subplot(2, 2, 1);
else
    subplot(1, 2, 1);
end
plot(s, Lambda_A, 'o-', s, Lambda_B, 's-');
xlabel('s_{tor}');
ylabel('coefficient [1]');
title('off-set');
legend('\Lambda_A (1/\surd\nu_*)', '\Lambda_B (1/\nu_*)');
grid on;

if has_sc
    subplot(2, 2, 2);
else
    subplot(1, 2, 2);
end
plot(s, Lambda_S, '^-');
xlabel('s_{tor}');
ylabel('coefficient [1]');
title('\Lambda_{S}');
grid on;

if has_sc
    subplot(2, 2, 3);
    plot(s, lambda_LC, 'o-');
    xlabel('s_{tor}');
    ylabel('coefficient [1]');
    title('\lambda^{LC}_{bB} (omnigenous Shaing-Callen)');
    grid on;

    subplot(2, 2, 4);
    plot(s, remainder, 's-');
    xlabel('s_{tor}');
    ylabel('coefficient [1]');
    title('remainder (non-omnigenous Shaing-Callen)');
    grid on;
end

set(gcf, 'Name', ['rabe output -- ' path]);
print('rabe_output.png', '-dpng', '-r150');
printf('saved rabe_output.png\n');
waitfor(gcf);
