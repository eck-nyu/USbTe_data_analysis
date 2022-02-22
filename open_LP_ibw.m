
%% 2022 Jan 10 
% First open Igor Pro from desktop
% Data --> Load Waves --> Load Igor Binary, select .ibw file, save as e.g. wave0
% Data --> Browse Waves, select the wave0 
% Wave Note contains analyzer info, photon energy, etc
% Data --> Save Waves --> Igor Text to save as .itx file
 
foldername = '/home/data/eck325/USbTe/2021dec_data/USbTe_1/LPolar/';
%% Import the LP data 
filenames =  {'LP7_98LH.itx'};
%{'LP11_98LH2.itx'};
% {'LP7_98LH.itx'}; ...
% 'LP8_98LV.itx';... 
% 'LP9_112LH.itx';... 
% 'LP10_112LV.itx';... 
% 'LP11_98LH2.itx';...
% 'LP12_125LV2.itx';... % LP12+ 2nd order harmonic % filename = 
% 'LP13_125LH.itx';... 
% 'LP14_132LH2.itx';...
% 'LP15_98LH2.itx';... % T=155K 
% 'LP16_112LV2_150K.itx';... % T=150K 
% 'LP17_112LV2_12K.itx';... % T=3~11K
% 'LP23_112LV_8K.itx'}; % T=8K

for fi = 1:numel(filenames)
    filename = filenames{fi};
[I, e, theta, polar] = load_itx_LP([foldername,filename]);
%% Take mean spectrum edc and set FL as steepest slope
spec = nanmean(I,3); 
edc = sum(spec,2);
% Interpolate to 3-meV step size 
eq = e(1): .003 : e(end); 
I = imresize( I, [numel(eq), numel(theta)]);
edcq = sum(nanmean(I,3),2);

% Check interpolation for artifacts (oscillations)
figure, plot(e, edc)
hold on, plot(eq, edcq,'r')
hold on, plot(0.5*(eq(2:end)+eq(1:end-1)), diff(edcq),'r');
hold on, plot(0.5*(e(2:end)+e(1:end-1)), diff(edc));

e = eq; edc = edcq; 

FL = e( find(diff(edc) == min(diff(edc))) );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FL_shift = 0.00;
FL = FL + FL_shift; % Manually adjust the FL 
disp(['File ',filename,', ','FL = ',num2str(FL),' eV']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Get constant-energy surface across theta(kx)-vs-polar(ky)

Ecut = -0.00; % Energy to see initial 2d surface
Ehw = 0.01; % Energy half-width to integrate 

thCutVal = 0.43;%0; % Theta val to see a sample spectral cut 
polCutVal = 7;%0; % Polar val to see a sample spectral cut

%%%%%%% Adjust the geometry angles/energies here %%%%%%%%%%%%%%%%%%%
azi0 = 0; tilt0 = -.5; pol0 = 0.8; hv = 98; WF = 4.4; EB = Ecut; V0 = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

geometry = [azi0, tilt0, pol0, hv, WF, V0];
[Esurf0, Kmap] = constantEsurf(I, e, theta, polar, FL, Ecut, Ehw, geometry);


% Plot constant-energy surface in both degrees, spec cuts at th,pol=0
figure,
subplot(223), imagesc(polar, theta, Esurf0), axis xy
hold on, plot(polCutVal*[1,1],[min(theta),max(theta)],'r--');
hold on, plot([min(polar),max(polar)],thCutVal*[1,1],'r--');
xlabel('polar (deg)'); ylabel('theta (deg)');
title(['Ecut = ',num2str(Ecut),' +/- ',num2str(Ehw),' eV']);

% Calculate index of theta,polar cuts and then plot the spectral cuts 
thCutIdx = round(interp1(theta, 1:numel(theta), thCutVal));
polCutIdx = round(interp1(polar, 1:numel(polar), polCutVal));

subplot(221), imagesc( polar, e-FL, reshape(I(:,thCutIdx,:), size(I,1),[]) ); axis xy;
hold on, plot([min(polar),max(polar)],(Ecut-Ehw)*[1,1],'r--');
hold on, plot([min(polar),max(polar)],(Ecut+Ehw)*[1,1],'r--');
title(['th = ',num2str(thCutVal),' deg']);

subplot(224), imagesc( e-FL, theta, I(:,:,polCutIdx)' ); axis xy; 
hold on, plot((Ecut-Ehw)*[1,1],[min(theta),max(theta)],'r--');
hold on, plot((Ecut+Ehw)*[1,1],[min(theta),max(theta)],'r--');
title(['pol = ',num2str(polCutVal),' deg']);

sgtitle(filename,'Interpreter','none'); colormap turbo;

% Plot 2d surfaces, at several sample energies, in kx-ky
Ecuts = [0, -0.1, -0.5, -1.0];
figure,
for E_i = 1:numel(Ecuts)
    Ecut = Ecuts(E_i);
    [Esurf, Kmap] = constantEsurf(I, e, theta, polar, FL, Ecut, Ehw, geometry);

    subplot(2,2,E_i);
    pcolor(Kmap(:,:,1), Kmap(:,:,2), Esurf), shading flat;
    hold on, plot([min(min(Kmap(:,:,1))),max(max(Kmap(:,:,1)))],0*[1,1],'w:');
    hold on, plot(0*[1,1],[min(min(Kmap(:,:,2))),max(max(Kmap(:,:,2)))],'w:');
    xlabel('Kx (invA)'); ylabel('Ky (invA)');
    title(['E = ',num2str(Ecut),' eV'])
end
sgtitle({filename},'Interpreter','none');
colormap turbo

%% Next step: Evaluate the gap at a given theta,polar (or kx,ky) coordinate 
thWidx = round( 0.5 / abs(theta(2)-theta(1)) ); % Width to integrate theta 
eStep = abs(mean(e(2:end)-e(1:end-1))); 
topX = 0.05; % Set energy peak as center-of-mass of top topX% of edc intensity 

polCoors = polar; % Polar values to get gap map 
thCoors = theta(thWidx+1:end-thWidx-1) ; % Theta values to get gap map 

thCoorsIdx = round(interp1( theta, 1:numel(theta), thCoors));
polCoorsIdx = round(interp1( polar, 1:numel(polar), polCoors));

gapMap = NaN*ones(numel(thCoorsIdx), numel(polCoorsIdx));
IdiffMap = NaN*ones( size(gapMap));

for th_i = 1:numel(thCoors)
    thIdx = thCoorsIdx(th_i);
    
    for pol_i = 1:numel(polCoors)
        polIdx = polCoorsIdx(pol_i); 
        
        spec = I( :, thIdx-thWidx : thIdx+thWidx, polIdx ); 
        [edc, e_centered] = symm_FL_edc(spec, e, FL );
        
        % Normalize edc for finding peak-com %
        edc = (edc-min(edc))/(max(edc)-min(edc));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        % Gap: find E of peak-com of symmetrized edc 
        % First find peak within 100meV of FL 
        findPkIdx = round(interp1( e_centered, 1:numel(e_centered), [-0.1, 0] )); 
        [pkHt, pkIdx] = max(edc(findPkIdx(1):findPkIdx(2))); 
        pkIdx = pkIdx + findPkIdx-1;
        
        % Find COM usin    sum(spec(:))g topX% of (normalized) edc intensity 
        topXstart = find( flipud(edc(1:pkIdx)) <= (1-topX)*pkHt, 1, 'first');
        topXstart = pkIdx - topXstart + 1;
        topXend = pkIdx-1 + find( edc(pkIdx:end) <= (1-topX)*pkHt, 1, 'first');        
        topXidx = topXstart : topXend; 
        thePkIdx = round( nansum( topXidx .* edc(topXidx)' )...
                            /nansum(edc(topXidx)) );
                        
        % Full gap is 2x the binding energy of the first peak in FL-symm edc  
        try
            theGap = 2 * abs( e_centered( thePkIdx )); % Using COM pk          
            gapMap( th_i, pol_i ) = theGap;
        end

        % NaN out values with inadequate statistics 
        bgndIdx = thePkIdx - round( 0.100 / eStep ); % Compare to bg int 100meV below band max
        try 
            Idiff = pkHt - edc(bgndIdx);       
            IdiffMap( th_i, pol_i ) = Idiff;
            
            if Idiff <= 0
                gapMap( th_i, pol_i ) = NaN;
            end
        end
    end
end

% Theta, polar grids for plotting with pcolor 
thetaMap = theta( repmat( reshape(thCoorsIdx, [],1), 1, numel(polCoors)));
polarMap = polar( repmat( reshape(polCoorsIdx, 1,[]), numel(thCoors), 1));

figure, pcolor(Kmap(thCoorsIdx,polCoorsIdx,1), Kmap(thCoorsIdx,polCoorsIdx,2), gapMap); shading interp;
% figure, pcolor(Kmap(thCoorsIdx,polCoorsIdx,1), Kmap(thCoorsIdx,polCoorsIdx,2), IdiffMap); shading interp;
% colormap(redblue); 
% xlabel('Polar (deg)'); ylabel('Theta (deg)');
xlabel('Kx (1/A)'); ylabel('Ky (1/A)'); 
colormap turbo; 
cb=colorbar(); ylabel(cb, '(eV)');
title({filename;'Gap map'},'Interpreter','none');
set(gca,'Color',[1,1,1]); % Set NaN values to white 


%% With a Kmap and gap map, interpolate into even kx,ky grid and start plotting kx/ky cuts

Kx = Kmap(thCoorsIdx, polCoorsIdx, 1); 
Ky = Kmap(thCoorsIdx, polCoorsIdx, 2);

KxLim = min(abs([min(Kx(:)),max(Kx(:))]));
Kxq = linspace(-KxLim, KxLim, numel(theta));
% Kxq = linspace(min(Kx(:)),max(Kx(:)), numel(theta));
Kyq = linspace(min(Ky(:)),max(Ky(:)), numel(polar));

GMq = griddata( reshape(Kx,1,[]), reshape(Ky,1,[]), reshape(gapMap,1,[]), Kxq, Kyq');
ESq = griddata( reshape(Kx,1,[]), reshape(Ky,1,[]), ...
        reshape(Esurf0(thCoorsIdx,polCoorsIdx),1,[]), Kxq, Kyq');
% There MUST be a way faster way to do this with a single-line function that can take unevenly-spaced 3d-data and interpolate. 
% For now since Kx,Ky are unevenly spaced, just use the griddata and do it for each energy 

Iq = zeros( numel(Kyq), numel(Kxq), numel(e) );
for e_i = 1:numel(e)
    
    Ie = reshape( I(e_i, thCoorsIdx, polCoorsIdx), numel(thCoorsIdx), numel(polCoorsIdx));
    Ieq = griddata( reshape(Kx,1,[]), reshape(Ky,1,[]), reshape(Ie,1,[]), Kxq, Kyq');    
    Iq(:,:,e_i) = Ieq;
end
        
Iqq = permute(Iq, [3,2,1]); 
Iqq = 0.5 * (Iqq + fliplr(Iqq));
%% Plot the gap map and couple Ky spectral cuts 
figure; 
subplot(4,3, [1,4]); 
pcolor(Kxq, Kyq, GMq); shading interp; xlabel('Kx'), ylabel('Ky') ; 
set(gca, 'color',[1,1,1]); colormap turbo; daspect([1,1,1]), 
ax1 = gca; hold(ax1, 'on');
title('Gap map'); 

subplot(4,3, [7,10]); 
pcolor(Kxq, Kyq, ESq); shading interp; xlabel('Kx'); ylabel('Ky');
colormap turbo; daspect([1,1,1]);
ax2 = gca; hold(ax2, 'on');
title('Fermi surface');

KyCuts = [0, -0.14, -0.30, -0.72, -1.0, -1.22, -1.38]; % Ky values to take spec cuts 
KyCutColors = hsv(numel(KyCuts)); 
subplotIdxs = [2,5,8,11,3,6,9,12]; % Subplot positions 

MDC_Es = [-0.1 : 0.05 : 0.05]; % Energies to plot MDCs 
EDC_Ks = [-1.0 : 0.2 : 1.0]; % Kx vals to plot EDCs 

MDC_Ehw = round( 0.01 / abs(e(2)-e(1))); % Energy range (idx) to +/- for MDCs
EDC_Khw = round( 0.05 / abs(Kxq(2)-Kxq(1))); % Kx range (idx) to +/- for EDCs

MDC_Eidxs = round(interp1(e-FL, 1:numel(e), MDC_Es));
EDC_Kidxs = round(interp1(Kxq, 1:numel(Kxq), EDC_Ks)); 

KyMDCs = zeros(numel(Kxq),numel(MDC_Es),numel(KyCuts)); % Array for storing MDCs
KyEDCs = zeros(numel(e),numel(EDC_Ks),numel(KyCuts)); % Array for storing EDCs

for Ky_i = 1:numel(KyCuts)
    KyCut = KyCuts(Ky_i);
    KyIdx = round( interp1( Kyq, 1:numel(Kyq), KyCut ));
    spec = Iqq(:,:,KyIdx);
    spec = spec / nansum(spec(:));

    % Get the MDCs and EDCs 
    for MDC_i = 1:numel(MDC_Es)
        Eidx = MDC_Eidxs(MDC_i);
        KyMDCs(:, MDC_i, Ky_i) = sum( spec( Eidx - MDC_Ehw : Eidx + MDC_Ehw, : ), 1);
    end
    for EDC_i = 1:numel(EDC_Ks)
        Kidx = EDC_Kidxs(EDC_i);
        KyEDCs(:,EDC_i,Ky_i) = sum( spec( :, Kidx - EDC_Khw : Kidx + EDC_Khw), 2);
    end
          
    % Plot the spectral cuts at each Ky value 
    subplot(4,3,subplotIdxs(Ky_i)); 
    imagesc(Kxq, e-FL, spec);
    hold on, plot([min(Kxq),max(Kxq)], 0*[1,1], 'w:');
    title(['Ky = ',num2str(KyCut)]);
    axis xy; colormap turbo; ylim([-0.5,0.1]);
    % Plot color-coded frame around subplot matching plot on gap map
    boxX = xlim; boxY = ylim; 
    hold on, rectangle('Position',[boxX(1),boxY(1),range(boxX),range(boxY)],'LineWidth',5,'EdgeColor', KyCutColors(Ky_i,:));
    
    % Plot the Ky cut on the gap map 
    plot(ax1, [min(Kxq),max(Kxq)], KyCut*[1,1], 'Color', KyCutColors(Ky_i,:))
    plot(ax2, [min(Kxq),max(Kxq)], KyCut*[1,1], 'Color', KyCutColors(Ky_i,:))
    
end
sgtitle(['FL shift = ',num2str(FL_shift),' eV']);
%%
figure;
t = tiledlayout(numel(KyCuts),2);
t.TileSpacing = 'Compact'; t.Padding = 'Compact'; 
for Ky_i = 1:numel(KyCuts)
   
    % Plot the spectral cuts at each Ky value 
    nexttile((Ky_i-1)*2+1)
    plot( repmat(Kxq',1,numel(MDC_Es)), linspace(0,1,numel(MDC_Es)) + KyMDCs(:,:,Ky_i) );
    legend(['Ky=',num2str(KyCuts(Ky_i))],'Location','northwest');
    
    nexttile((Ky_i-1)*2+2)
    plot( linspace(0,1,numel(EDC_Ks)) + KyEDCs(:,:,Ky_i), repmat(e'-FL,1,numel(EDC_Ks)) );
end
end


