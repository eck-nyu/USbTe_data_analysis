
load('spec3d_LHLV.mat');
azi0 = 0; tilt0 = 0.45; pol0 = 0; WF = 4.4; V0 = 13; 
geometry_angles = [azi0,tilt0,pol0,WF,V0];

Ecuts = [0:-.5:-1];
c = 1.0*9.063; % Angstrom
a = 1.0*4.321; % Angstrom

% spec3d = spec3dLH + spec3dLV;

% Plot the kx-kz dispersion at several energies with grids 
fig = figure('Renderer','Painters'); 
for Ecut_i = 1:numel(Ecuts)
    Ecut = Ecuts(Ecut_i);
    
    [kxMap,kzMap,specCut] = hv2kz( spec3d, Ecut, E, th, T, geometry_angles, a, c);
    kxMap = kxMap(:, 1:7:end);
    kzMap = kzMap(:, 1:7:end); 
    specCut = imresize(specCut, [size(kxMap,1), size(kxMap,2)]); 
    hold on, surf( kxMap, kzMap, Ecut*ones(size(kxMap)), specCut ), shading interp
    
     % Plot the reference hv values 
    hold on, text(kxMap(1:3:end,end), kzMap(1:3:end,end), Ecut*ones(size(kxMap(1:3:end,end))), strsplit(num2str(round(T(1:3:end)))))
    
    kzLim = ylim; 
    for nkz_i = ceil(kzLim(1)):floor(kzLim(2))
        if rem(nkz_i,2) == 0
            hold on, plot3([min(kxMap(:)),max(kxMap(:))], nkz_i*[1,1], Ecut*[1,1],'r:');
        else
            hold on, plot3([min(kxMap(:)),max(kxMap(:))], nkz_i*[1,1], Ecut*[1,1],'r-');
        end
    end
    kxLim = xlim; 
    for nkx_i = -2:2
        if rem(nkx_i,2)==0
            hold on, plot3(nkx_i*[1,1],[min(kzMap(:)),max(kzMap(:))], Ecut*[1,1], 'r:');
        else 
            hold on, plot3(nkx_i*[1,1],[min(kzMap(:)),max(kzMap(:))], Ecut*[1,1], 'r-');
        end  
    end
    title({['lattice a,c: [',num2str(round(a,2)),', ',num2str(round(c,2)),'] A'];['V0: ',num2str(geometry_angles(5)),' eV']})
    

end
colormap turbo
xlabel('kx/(\pi/a)'); ylabel('kz/(\pi/c)'); zlabel('E (eV)');
view(10,30), grid on;
caxis([0.0005,0.0065])
%%
figure,
hv_is = 1:5:numel(T);
for hv_i = 1:numel(hv_is)
    hv_ii = hv_is(hv_i);
    hv = T(hv_ii);
    EB = 0;
    kx = theta2kx(th,azi0,tilt0,pol0,hv,WF,EB,V0);
    subplot(numel(hv_is),1,hv_i);
    imagesc(kx/(pi/a), E, spec3d(:,:,hv_ii)); axis xy, 
    title(['hv=',num2str(hv),'eV']);
    ylim([-.5,.05]);
end
colormap turbo

function [kxMap,kzMap,EcutSpec] = hv2kz( spec3d, Ecut, Evec, thetaVec, ...
                                        hvVec, geometry_angles, a,c ) 
% geometry_angles vector 1.azi0, 2.tilt0, 3.pol0, 4.WF, 5.V0 
plot_figure = 0; 

azi0 = geometry_angles(1); 
tilt0 = geometry_angles(2); 
pol0 = geometry_angles(3); 
WF = geometry_angles(4); 
V0 = geometry_angles(5); 
EB = -Ecut; % Approximation: Assume EB is 0 for all energies

% Eidx = round(interp1( Evec, 1:numel(Evec), Ecut));
Ewidx = round( .02 / abs(Evec(2)-Evec(1)) );

% EcutSpec = sum(spec3d( Eidx-Ewidx:Eidx+Ewidx, :, :), 1);
% EcutSpec = permute(EcutSpec, [3,2,1]); % Permute from E,th,hv to hv,th(,E) 

EcutSpec = zeros(numel(hvVec), numel(thetaVec));
kxMap = zeros(numel(hvVec), numel(thetaVec));
kzMap = zeros(numel(hvVec), numel(thetaVec)); 
for hv_i = 1:numel(hvVec)
    hv = hvVec(hv_i);
    kMap = theta2kMap( thetaVec, azi0, tilt0, pol0, hv, WF, EB, V0);
    
    kxMap(hv_i,:) = kMap(1,:);
    kzMap(hv_i,:) = kMap(3,:);
    
    % Calculate FL for each hv
    spec = spec3d(:,:,hv_i); 
    edc = sum(spec, 2);
    FL_range = find(Evec>-.05,1,'first') : find(Evec>0.05,1,'first'); % Assume FL not off more than -1,.05 eV 
    slope = diff(edc);
    FL_idx = FL_range(1)-1+find( slope(FL_range) == min(slope(FL_range)),1,'first');
    FL = Evec( FL_idx );
    theEvec = Evec - FL; 
    
    Eidx = round(interp1( theEvec, 1:numel(theEvec), Ecut ));
    
    EcutSpec_ = sum(spec3d( Eidx-Ewidx:Eidx+Ewidx, :,hv_i),1);
    EcutSpec_ = permute(EcutSpec_,[3,2,1]); 
    EcutSpec(hv_i,:) = EcutSpec_ / sum(EcutSpec_);
end

% Plot the kx-kz dispersion
kxMap = kxMap / (pi/a); % Rescale kx,kz in units of kz,kc 
kzMap = kzMap / (pi/c);

if plot_figure == 1
    figure,

    pcolor(kxMap,kzMap,EcutSpec), shading interp, colormap turbo
    daspect([a/pi, c/pi,1])

    % Plot the reference hv values 
    hold on, text(kxMap(:,end),kzMap(:,end),strsplit(num2str(round(hvVec))))

    % Plot horizontal lines for integer-kz 
    if ~isempty(c)
        kzLim = ylim; 
        for nkz_i = floor(kzLim(1)) : floor(kzLim(2))
            if rem(nkz_i,2) == 0
                hold on, plot([min(kxMap(:)),max(kxMap(:))], nkz_i*[1,1],'r:');
            else
                hold on, plot([min(kxMap(:)),max(kxMap(:))], nkz_i*[1,1],'r-');
            end
        end
        kxLim = xlim; 
        for nkx_i = floor(kxLim(1)) : floor(kxLim(2))
            if rem(nkx_i,2)==0
                hold on, plot(nkx_i*[1,1],[min(kzMap(:)),max(kzMap(:))], 'r:');
            else 
                hold on, plot(nkx_i*[1,1],[min(kzMap(:)),max(kzMap(:))], 'r-');
            end  
        end
    end
    title({['lattice a,c: [',num2str(round(a,2)),', ',num2str(round(c,2)),'] A'];['V0: ',num2str(V0),' eV'];['E= ',num2str(Ecut),' eV']})
    caxis(caxRange(EcutSpec, .1,.95))
    xlabel('kx /(\pi/a)'); 
    ylabel('kz /(\pi/c)');
end
% pcolor(kxMap,kzMap,EcutSpec), shading interp, colormap turbo
% daspect([1,1,1]);
% 
% % Plot the reference hv values 
% hold on, text(kxMap(:,end),kzMap(:,end),strsplit(num2str(round(hvVec))))
% 
% % Plot horizontal lines for integer-kz 
% if ~isempty(c)
%     kzLim = ylim % range of kz
%     nkzLim = kzLim / (pi/c) % integer units of kz
%     nkzTicks = floor(nkzLim(1)):1:floor(nkzLim(2));
%     for nkz_i = 1:numel(nkzTicks)
%         nkzTick = nkzTicks(nkz_i) * (pi/c);
%         if rem(nkzTicks(nkz_i),2) == 0
%             hold on, plot([min(kxMap(:)),max(kxMap(:))], nkzTick*[1,1],'r--');
%         else
%             hold on, plot([min(kxMap(:)),max(kxMap(:))], nkzTick*[1,1],'r:');
%         end
%     end
%     nkxLim = xlim / (pi/a);
%     nkxTicks = floor(nkxLim(1)) : 1 : floor(nkxLim(2));
%     for nkx_i = 1:numel(nkxTicks)
%         nkxTick = nkxTicks(nkx_i) * (pi/c);
%         if rem(nkxTicks(nkz_i),2)==0
%             hold on, plot(nkxTick*[1,1],[min(kzMap(:)),max(kzMap(:))], 'r--');
%         else 
%             hold on, plot(nkxTick*[1,1],[min(kzMap(:)),max(kzMap(:))], 'r:');
%         end  
%     end
% end

end

% 
% kMap3d = zeros( numel(Evec), numel(thetaVec), numel(hvVec), 3 );
% for hv_i = 1:numel(hvVec)
%     hv = hvVec(hv_i);
%     kMap = theta2kMap( thetaVec, azi0, tilt0, pol0, hv, WF, EB, V0);
%     
%     kxMapE = repmat(kMap(1,:), numel(Evec), 1); % Repeat kMap for all energies
%     kyMapE = repmat(kMap(2,:), numel(Evec), 1); 
%     kzMapE = repmat(kMap(3,:), numel(Evec), 1);
%     
%     kMap3d(:,:,hv_i,1) = kxMapE; 
%     kMap3d(:,:,hv_i,2) = kyMapE; 
%     kMap3d(:,:,hv_i,3) = kzMapE; 
% end