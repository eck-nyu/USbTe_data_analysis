%% Demo script for opening 2021-Dec ALS beamtime data from JD
% Theta is analyzer slit axis, polar is sample surface tilt axis 

%% Open single-cut (energy-vs-theta) files 
foldername = '/home/data/eck325/USbTe/2021dec_data/USbTe_1/itx/';
filename = 'ust_046_S001.itx'; % T=134K, hv=112eV(LV), pol=0deg
filename = 'ust_039_S036.itx'; % T=155K, hv=98eV(LH), pol=7deg
filename = 'ust_082.itx'; % Au ref, hv=98LH, T~8K

[spec, energy, theta] = load_itx([foldername,filename]); % function loads single-cut itx file into 2d image with e,theta axes

% Plot the spectrum 
figure, imagesc(theta, energy, spec); % Omit "axis xy" here since these particular scans have energy in units of binding energy
xlabel('theta (deg)'); ylabel('binding energy (eV)')
%% Open LivePolar scan (energy-vs-theta-vs-polar) files
filename = 'LP16_112LV2_150K.itx'; % T=150K, hv=112eV(LV), pol=-14:28deg
filename = 'LP15_98LH2.itx'; % T=155K, hv=98eV(LH), pol=-14:28deg

[specs, e, theta, polar] = load_itx_LP([foldername,filename]); % function loads the LP file into a 3d image array with e,theta,polar axes

% Estimate FL for now using edc slope minimum 
edc = sum(sum(specs, 3),2);
FL = e( find( diff(edc) == min(diff(edc))) );

figure, plot(e, edc); xlabel('energy (eV)'); ylabel('all-int edc'); 
yyaxis right, plot(0.5*(e(2:end)+e(1:end-1)), diff(edc)), ylabel('edc slope');

Ecut = -0.0; % Energy to see constant-E surface (not binding energy)
Ehw = 0.05; % Sum Ecut +/- Ehw
cut_theta_val = 0; % See a spec cut through theta value
cut_polar_val = 7; % See a spec cut through polar value 

hv = 98; %%%%%%%%%%%%%%%%%%% UPDATE THIS VALUE
azi0=0; tilt0=0; WF=4.4; innerPot=10; % Quick approx. placeholders for exp setup values
geometry = [azi0, tilt0, hv, WF, innerPot];

[Esurf, Kmap] = constantEsurf(I, e, theta, polar, FL, Ecut, Ehw, geometry); % Get the 2d image array of constant-energy surface

 
figure,
% Plot constant-energy surface at specified Ecut, Ehw using estimate FL 
subplot(223), imagesc(polar, theta, Esurf), axis xy
hold on, plot(cut_polar_val*[1,1],[min(theta),max(theta)],'r--');
hold on, plot([min(polar),max(polar)],cut_theta_val*[1,1],'r--');
xlabel('polar (deg)'); ylabel('theta (deg)');

% Plot spec cuts 
th0idx = round(interp1(theta, 1:numel(theta), cut_theta_val)); % find theta index
pol0idx = round(interp1(polar, 1:numel(polar), cut_polar_val)); % find polar index

subplot(221), imagesc( polar, e, reshape(specs(:,th0idx,:), size(specs,1),[]) ); axis xy;
hold on, plot( [min(polar),max(polar)], (FL+Ecut+Ehw)*[1,1], 'r--');
hold on, plot( [min(polar),max(polar)], (FL+Ecut-Ehw)*[1,1], 'r--');
xlabel('polar (deg)'); ylabel('energy (eV)');

subplot(224), imagesc( e, theta, specs(:,:,pol0idx)' ); axis xy; 
hold on, plot((FL+Ecut+Ehw)*[1,1], [min(theta),max(theta)], 'r--');
hold on, plot((FL+Ecut-Ehw)*[1,1], [min(theta),max(theta)], 'r--');
ylabel('theta (deg)'); xlabel('energy (eV)');
colormap turbo