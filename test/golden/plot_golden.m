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
Lambda_bl     = ncread(path, 'Lambda_bl');
Lambda_lm     = ncread(path, 'Lambda_lm');
Lambda_finite = ncread(path, 'Lambda_finite');

info = ncinfo(path);
varnames = {info.Variables.Name};
has_sc = any(strcmp(varnames, 'lambda_SC_bB'));

if has_sc
    lambda_SC = ncread(path, 'lambda_SC_bB');
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
plot(s, Lambda_bl, 'o-', s, Lambda_lm, 's-');
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
plot(s, Lambda_finite, '^-');
xlabel('s_{tor}');
ylabel('coefficient [1]');
title('\Lambda_{HGM}');
grid on;

if has_sc
    subplot(2, 2, 3);
    plot(s, lambda_SC, 'o-');
    xlabel('s_{tor}');
    ylabel('coefficient [1]');
    title('\lambda^{SC}_{bB} (omnigenous Shaing-Callen)');
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
