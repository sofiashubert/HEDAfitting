function [] = HEDA_MultiFit()

fileName = importdata('test.txt');

% Get the experimental extinction
omega = fileName(:,1);
ext_data = fileName(:,2:end);

% sort data so that it's increasing
if omega(1)>omega(2)
     omega = flip(omega);
     ext_data = flip(ext_data);
end


% Constants
m_0 = 9.109e-31; % mass of an electron (kg)
m_eff_init = 0.38; % initial effective mass
eps_0 = 8.854e-12; % vacuum permittivity (F/m)
mu_0 = 1.257e-6; % vacuum permeability (H/m)
e = 1.602e-19; % electron charge (C)
hbar = 1.055e-34; % Plank's constant (J s)
c = 1/(sqrt(eps_0*mu_0)); % speed of light (m/s)
C_nonparabolic = 0.5;  % Nonparabolicity coefficient (CdO = 0.5)
hbar_eV = 6.582119569e-16;   % Plank's constant units (eV s)
A_surf = 0.25; % proportionality factor of surface damping (CdO = 0.25, ITO = 0.75)

%%%%%%%%%%%%%%  Experimental parameters: change with each titration  %%%%%%%%%%%%%%
CoCp2_vol = [0  0	5	10	15	20	25	30	35	40	45	50]; % Name of spectra
a_mean = 1/2*[16.1	16.1	16.1	16.1	16.1	16.1	16.1	16.1	16.1	16.1	16.1	16.1]*10^-9; % mean nanoparticle radius (m)
a_dev = 1/2*[1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2	1.2]*10^-9; % radius standard deviation (m)

eps_inf = 5.5*eps_0; % high-freq permittivity (F/m) (CdO = 5.5)
eps_f = 1.4515^2*eps_0; % solvent permittivity (F/m) (THF+tol = 1.4525)
path_length = 0.001; % path length (m)
eta = [5.19653E-05	4.45E-05	4.39E-05	4.33E-05	4.27E-05	4.21E-05	4.16E-05	4.10E-05	4.05E-05	4.00E-05	3.95E-05	3.90E-05]; % NC volume fraction from ICP

% Initializations
m_e_init = m_eff_init*m_0; % initialize the rescaling of m_e with Tauc m_e
omega_LSPR = zeros(size(CoCp2_vol));
ext_LSPR = zeros(size(CoCp2_vol));
w_LSPR = zeros(size(CoCp2_vol));
ext_HEDA = zeros(length(omega), length(CoCp2_vol));
gamma_mean_vals = zeros(size(CoCp2_vol));
maxY_vals = zeros(size(CoCp2_vol));
maxX_vals = zeros(size(CoCp2_vol));
ne_vals = zeros(size(CoCp2_vol));
ne_std_vals = zeros(size(CoCp2_vol));
fe_vals = zeros(size(CoCp2_vol));
ell_vals = zeros(size(CoCp2_vol));
m_e_vals = zeros(size(CoCp2_vol));
m_e_vals(:,1) = m_eff_init;


close all
h_all = figure;
hold on
colors = parula(length(CoCp2_vol));

    function ext = Drude_fun(omega, omega_p_f, gamma_tilde, eta, path_length)

    % INPUTS
    % omega = (Nk-by-1) frequency (cm^-1)
    % omega_p_f = (scalar) plasma frequency in fluid (cm^-1)
    % gamma_tilde = (scalar) damping (relative to omega_p_f)
    % eta = (scalar) particle volume fraction
    % path_length = (scalar) path length
    %
    % 
    % OUTPUTS
    % ext = (Nk-by-1) extinction


    % Permittivity
    eps_p_tilde = eps_inf/eps_f - 1./((omega/omega_p_f).^2+1i*(omega/omega_p_f).*gamma_tilde);

    % Polarizability
    alpha = (eps_p_tilde-1)./(eps_p_tilde+2); % polarizability 

    % Extinction cross-section per NP volume (1/m)
    omega =  2*pi*c*100*omega; % convert cm^-1 to s^-1
    ext_coeff = 3*sqrt(eps_f*mu_0).*omega.*imag(alpha);

    % Total extinction -log(I_t/I_0)
    ext = eta*path_length.*ext_coeff/log(10);

end

    function ext = ext_fun(omega, n_mean, n_dev, ell, eta_c, a_mean, a_dev, eta, path_length,i)

    % INPUTS
    % omega = frequency (a^-1)
    % n_mean = mean carrier conc. (a^-3)
    % n_dev = std. dev. carrier conc. (a^-3)
    % ell = mean free path (a)
    % eta_c = core volume fraction
    % a_mean = particle radius (m)
    % a_dev = std. dev. particle radius (m)
    % eta = (scalar) particle volume fraction
    % path_length = (scalar) path length (m)
    % m_eff_init = (scalar) Tauc based m_e to use
    % i = spectrum iteration number
    %
    % OUTPUTS
    % ext = extinction
    % 

    % Convert the inputs to dimensional form
    omega = omega/a_mean;
    n_mean = n_mean/a_mean^3;
    n_dev = n_dev/a_mean^3;
    ell = ell*a_mean;
    
    % Create grid of carrier conc. and radii
    n_lin = linspace(n_mean-4*n_dev, n_mean+4*n_dev, 101); % carrier conc. out several standard devs from the mean
    a_lin = linspace(a_mean-4*a_dev, a_mean+4*a_dev, 101); % radius out several standard devs from the mean
    
    % if fitting the as-synthesized spectrum, use the initial effective
    % mass from Tauc analysis
    % create complementary m_e matrix calculated from n_lin matrix   
    if i == 1
        m_e_lin = linspace(m_e_init,m_e_init,101);    
    else
        coeff_m_e = 0.16*m_0*sqrt(1+(2*C_nonparabolic*((hbar_eV^2)/m_0)*(3*pi^2)^(2/3)));  
        m_e_lin = linspace(coeff_m_e*((n_mean/10^27)^(1/3)),coeff_m_e*((n_mean/10^27)^(1/3)),101);
    end

    %[n, m_e, a] = ndgrid(n_lin, m_e_lin, a_lin); % changed this to 3D grid: carrier conc., effective mass , radius % probably the problem
    [n, m_e] = ndgrid(n_lin, m_e_lin);
    [n, a] = ndgrid(n_lin, a_lin);
 

    % Probability density function for the carrier conc. and radii
    P = 1/(2*pi*n_dev*a_dev)*exp(-1/2*((n-n_mean).^2/n_dev^2 + (a-a_mean).^2/a_dev^2));
    
    % freq. gets the third dimension
    omega = permute(omega, [3, 2, 1]);
   
    % Core dielectric function
    omega = 2*pi*c*omega; % convert m^-1 to s^-1
    omega_p = sqrt(e^2*n./(eps_0*m_e)); % core plasma freq. (s^-1)
    gamma = (3*pi^2)^(1/3)*hbar./m_e.*n.^(1/3).*((A_surf./(a.*eta_c.^(1/3)))+1./ell); % damping coefficient (s^-1) 
    eps_c = eps_inf-eps_0*omega_p.^2./(omega.^2+1i*omega.*gamma); % core dielectric function relative to fluid
   
    % Shell dielectric function
    eps_s = eps_inf;

    % Core-shell dielectric function relative to fluid
    eps_p = eps_s.*((eps_c+2*eps_s)+2*eta_c.*(eps_c-eps_s))./((eps_c+2*eps_s)-eta_c.*(eps_c-eps_s));

    % Extinction
    alpha = (eps_p/eps_f-1)./(eps_p/eps_f+2); % polarizability 
    ext_coeff = 4*pi*a.^3.*sqrt(eps_f*mu_0).*omega.*imag(alpha); % extinction cross-section (m^2)
    
    % Ensemble-averaged values
    ext_coeff_ave = squeeze(trapz(a_lin, trapz(n_lin, ext_coeff.*P, 1), 2)); % extinction cross-section (m^2)
    V_ave = trapz(a_lin, trapz(n_lin, 4*pi*a.^3/3.*P)); % NP volume (m^3)
    ext = eta*path_length.*ext_coeff_ave/(log(10)*V_ave); % extinction
    
end



% Loop through doping amounts
for i = 1:(length(CoCp2_vol))

    % Current extinction
    ext = ext_data(:,i);
    
    % Get the full width half max
    [~,~,~,p] = findpeaks(ext, omega); % find the prominences of all peaks
    [ext_LSPR(i), omega_LSPR(i), w_LSPR(i), ~] = findpeaks(ext, omega, 'MinPeakProminence', max(p), 'WidthReference', 'halfheight'); % work with only the most prominent peak
    
    % Use this to get good initial guesses for a simple Drude fit
    omega_p_f_guess = omega_LSPR(i)*sqrt(eps_inf/eps_f+2); % (cm^-1)
    gamma_tilde_guess = eta(i)*path_length/(ext_LSPR(i)*log(10))*9*sqrt(eps_f*mu_0)*(2*pi*c*100*omega_p_f_guess)/(eps_inf/eps_f+2)^2; % (relative to omega_p_f)

    % Fit the data only for frequecies in a range of w_LSPR centered
    % around omega_LSPR
    range = 1.2;
    omega_lim = [omega_LSPR(i) - range*w_LSPR(i)/2, omega_LSPR(i) + range*w_LSPR(i)/2];
    ind = (omega >= omega_lim(1)) & (omega <= omega_lim(2));
    
    omega_fit = omega(ind);
    ext_fit = ext(ind);
    
    % Start with a simple Drude fit to get a good initial guess
    ft = fittype(@(omega_p_f, gamma_tilde, eta, omega) Drude_fun(omega, omega_p_f, gamma_tilde, eta, path_length), 'independent', 'omega', 'dependent', 'ext');
    opts = fitoptions(ft);
    opts = fitoptions(opts, 'Lower', [0, 0, 0], 'Upper', [Inf, Inf, 1], 'StartPoint', [omega_p_f_guess, gamma_tilde_guess, eta(i)]);
    f_Drude = fit(omega_fit, ext_fit, ft, opts);
    coeffs = coeffvalues(f_Drude);

    % Get guesses for the HEDA parameters from the Drude fit
    n_mean_guess = (2*pi*c*100*coeffs(1))^2*eps_f*m_e_init/e^2; % (m^-3)
    gamma_mean_guess = coeffs(2)*(2*pi*c*100*coeffs(1)); % (m^-3)
    ell_guess = (m_e_init*gamma_mean_guess/((3*pi^2)^(1/3)*hbar*n_mean_guess^(1/3)) - 3/(4*a_mean(i)))^-1; % (m)
    
    % Convert to dimensionless form
    n_mean_guess = n_mean_guess*a_mean(i)^3; % (a^-3)
    n_dev_guess = 0.1*n_mean_guess; % (a^-3)    
    eta_c_guess = 0.5;

    % Choose upper bound for ell and recompute guess if it is out of bounds
    % 100 upper bound for ICO, check bulk mfp for other materials
    ell_bound = 100*10^-9/a_mean(i);
    if (ell_guess >= ell_bound) || (ell_guess < 0)
        ell_guess = 0.9*ell_bound;
    end
    
    % Fit the HEDA model to the data
    ft = fittype(@(n_mean, n_dev, ell, eta_c, omega) ext_fun(omega, n_mean, n_dev, ell, eta_c, a_mean(i), a_dev(i), eta(i), path_length,i), 'independent', 'omega', 'dependent', 'ext');
    opts = fitoptions(ft);
    opts = fitoptions(opts, 'Lower', [0, 0, 0, 0], 'Upper', [Inf, Inf, ell_bound, 1], 'StartPoint', [n_mean_guess, n_dev_guess, ell_guess, eta_c_guess]);
    f_HEDA = fit(omega_fit*100*a_mean(i), ext_fit, ft, opts);
    coeffs = coeffvalues(f_HEDA);

    % Calculate the avg. m_e used based on the output mean n_e
    if i ==1 
        m_e_final = m_eff_init;
    else
        m_e_final = 0.16*sqrt(1+(2*C_nonparabolic*(((hbar_eV)^2)/(m_0)*(3*pi^2*((coeffs(1)/a_mean(i)^3)/10^27))^(2/3))));
    end
    
    % Compute some things from the fit and append to array with all values
    omega_p = sqrt(e^2*coeffs(1)/a_mean(i)^3/(eps_0*(m_e_final*m_0)))/(2*pi*100*c); % plasma freq (cm^-1)
    omega_p_f = omega_p*sqrt(eps_0/eps_f); % plasma freq in fluid (cm^-1)
    
    gamma_mean = (3*pi^2)^(1/3)*hbar/(m_e_final*m_0)*(coeffs(1)/a_mean(i)^3).^(1/3).*(A_surf./(a_mean(i).*coeffs(4).^(1/3))+1./(coeffs(3)*a_mean(i)))/(2*pi*100*c); % damping (cm^-1)
    gamma_mean_vals(:,i) = gamma_mean;

    ne_vals(:,i) = coeffs(1)/a_mean(i)^3;
    ne_std_vals(:,i) = coeffs(2)/a_mean(i)^3;
    fe_vals(:,i) = coeffs(4);
    ell_vals(:,i) = coeffs(3)*a_mean(i);
    m_e_vals(:,i) = m_e_final;
 
    %Calculate max X and Y values
    maxY_vals(:,i) = max(ext_fit(:,1));
    indexY= find(ext_fit(:,1) == maxY_vals(1,i));
    maxX = omega_fit(indexY);
    maxX_vals(:,i) = maxX;


    % Print results -- useful if you're doing a single fit but writing to
    % excel file is better for titration data
    fprintf('HEDA fit:\n')
    fprintf('CoCp*2 vol = %.2f, a = %.2e nm\n', CoCp2_vol(i), a_mean(i))
    fprintf('n_mean = %.4f = %.3e m^-1/3\n', coeffs(1), coeffs(1)/a_mean(i)^3)
    fprintf('n_dev = %.4f = %.3e m^-1/3\n', coeffs(2), coeffs(2)/a_mean(i)^3)
    fprintf('omega_p = %.0f cm^-1\n', omega_p)
    fprintf('omega_p_f = %.0f cm^-1\n', omega_p_f)
    fprintf('gamma/omega_p_f = %.3f, gamma = %.0f cm^-1\n', gamma_mean/omega_p_f, gamma_mean)
    fprintf('ell = %.3f = %.3g m\n', coeffs(3), coeffs(3)*a_mean(i))
    fprintf('eta_c = %.3f\n\n', coeffs(4))


    %Write HEDA results to an excel file
    results_ColumnNames = {'n_e','n_e_std','fe','mfp','gamma','m_eff','maxFreq','maxExt'};
    T_FitResults = table(transpose(ne_vals/10^26),transpose(ne_std_vals/10^26),transpose(fe_vals),transpose(ell_vals*10^9),transpose(gamma_mean_vals),transpose(m_e_vals),transpose(maxX_vals),transpose(maxY_vals),'VariableNames', results_ColumnNames);
    
    %%%%%%%%%%%%%%  change with each titration  %%%%%%%%%%%%%%
    writetable(T_FitResults, 'test_HEDA_FitResults.xlsx');

     % Plot all experimental and fit data
    figure(h_all)
    plot(omega, ext, 'Color', colors(i,:), 'LineWidth', 3)
    plot(omega, f_HEDA(omega*100*a_mean(i)), ':',  'Color', colors(i,:), 'LineWidth', 3)

    xlabel('\omega (cm^{-1})')
    ylabel('extinction')
    xlim([1500, 7000])
    set(gca,'FontUnits','normalized','FontSize',0.05,...
        'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
    set(gca, 'XDir', 'reverse')
    
    % Store fits
    ext_HEDA(:,i) = f_HEDA(omega*100*a_mean(i));


figure(h_all)
xlabel('\omega (cm^{-1})')
ylabel('extinction')
xlim([4000, 7000])
set(gca,'FontUnits','normalized','FontSize',0.05,...
    'FontWeight','bold','LineWidth',1,'PlotBoxAspectRatio',[1,1,1])
set(gca, 'XDir', 'reverse')

% Save fits to an Excel file
%column_name = {'freq.', '1 (data)', '2 (data)', '3 (data)', '4 (data)', '5 (data)', '6 (data)', '1 (fit)', '2 (fit)', '3 (fit)', '4 (fit)', '5 (fit)', '6 (fit)'};
column_name = {'Wavenumber', 'NC_AS','NC_acid','Aliq1','Aliq2','Aliq3','Aliq4','Aliq5','Aliq6','Aliq7','Aliq8','Aliq9','Aliq10'};
%T_FitExt = table(omega, ext_data(:,1), ext_data(:,2), ext_data(:,3), ext_data(:,4), ext_data(:,5), ext_data(:,6), ext_HEDA(:,1), ext_HEDA(:,2), ext_HEDA(:,3), ext_HEDA(:,4), ext_HEDA(:,5), ext_HEDA(:,6), 'VariableNames', column_name);
T_FitExt = table(omega,ext_HEDA(:,1), ext_HEDA(:,2), ext_HEDA(:,3), ext_HEDA(:,4), ext_HEDA(:,5), ext_HEDA(:,6),ext_HEDA(:,7), ext_HEDA(:,8), ext_HEDA(:,9), ext_HEDA(:,10), ext_HEDA(:,11), ext_HEDA(:,12), 'VariableNames', column_name);

%%%%%%%%%%%%%%  change with each titration  %%%%%%%%%%%%%%
writetable(T_FitExt, 'test_HEDA_FitExt.xlsx')
end


end