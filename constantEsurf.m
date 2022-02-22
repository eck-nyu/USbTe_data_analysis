
function [Esurf, Kmap] = constantEsurf( I, e, theta, polar, FL, Ecut, Ehw, geometry )
% Update 2022 01 21: Added element to geometry array, pol0 to be scalar
% offset added to polar vector, account for extra polar angle between manipulator and
% sample surface

% Take 3d image data (LP scan), input an energy value for constant-energy
% cuts to plot a 2d constant energy surface 

% INPUTS
% I is 3d image data, dimensions 1=energy (e), 2=theta (analyzer slit
% angle), 3=polar (manipulator surface tilt). 
% FL is energy value of fermi level
% Ecut in units of negative binding energy (E-FL)
% Ehw is half-width of range to sum Ecut +/- Ehw
% geometry is experimental setup angles, energies 

% OUTPUTS 
% Esurf 2d array constant energy surface, dimensions 1=theta, 2=polar 
% Kmap(:,:,1) is 2d array of corresponding kx values
% Kmap(:,:,2) is 2d array of corersponding ky values

azi0 = geometry(1); % azimuth 
tilt0 = geometry(2); % tilt 
pol0 = geometry(3); % polar offset 
hv = geometry(4); % photon energy 
WF = geometry(5); % Work function
innerPot = geometry(6); % Inner potential 

EB = -Ecut; % binding energy 

Estep = abs(mean(e(2:end)-e(1:end-1)));
% Convert Ecut, Ehw to index
EcutIdx = round(interp1( e-FL, 1:numel(e), Ecut));
EhwIdx = round(Ehw / Estep);

Esurf = zeros(numel(theta), numel(polar));
Kmap = zeros(numel(theta), numel(polar), 3);
for i = 1:numel(polar)
    spec = I(:,:,i);
    
    mdc = sum(spec(EcutIdx - EhwIdx : EcutIdx + EhwIdx, :), 1);   
    Esurf(:,i) = mdc;
    
    kMap = theta2kMap( theta, azi0, tilt0, polar(i)+pol0, hv, WF, EB, innerPot);
    Kmap(:,i,1) = kMap(1,:);
    Kmap(:,i,2) = kMap(2,:);
end


end