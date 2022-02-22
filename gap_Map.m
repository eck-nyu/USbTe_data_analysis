
%% 2022 Jan 10 
% First open Igor Pro from desktop
% Data --> Load Waves --> Load Igor Binary, select .ibw file, save as e.g. wave0
% Data --> Browse Waves, select the wave0 
% Wave Note contains analyzer info, photon energy, etc
% Data --> Save Waves --> Igor Text to save as .itx file
 
foldername = '/home/data/eck325/USbTe/2021dec_data/USbTe_1/LPolar/';
%% Import the LP data 
filenames = {'LP7_98LH.itx'}; %{'LP15_98LH2.itx'};%{'LP11_98LH2.itx'};
                combine_LH_and_LV = 1; 
                % To combine LH,LV scans, run LV first, assign ILV = I.
                % Then set this to 1 and run with LH data 
                
%%%%%%% Adjust the geometry angles/energies here %%%%%%%%%%%%%%%%%%%
% list of filename IDs and azi0,tilt0,pol0,FL_shift 
% Build up the parameter table as you go. Play with geometry values and
% look at plots, store them as you find them. 
params = table('Size',[1,6],'VariableNames',["ID","hv","azi0","tilt0","pol0","FLshift"],'VariableTypes',["string","double","double","double","double","double"]);
params(1,:)     = {'LP7',98,  -0.4, -0.45, 1.0, 0};
params(end+1,:) = {'LP8',98,  -0.4, -0.65, 1.0, 0};
params(end+1,:) = {'LP9',112,  0.0, -0.50, 0.5, 0.003};
params(end+1,:) = {'LP10',112, 0.0, -0.50, 0.5, 0.003};
params(end+1,:) = {'LP15',98,  1.0,  0.00, 0.5, 0};
params(end+1,:) = {'LP17',112, 1.0,  0.00, 0.5, 0};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KyCuts = [0.76, 0.255, 0, -0.68, -1.45]; 
 

for fi = 1:numel(filenames)
filename = filenames{fi};
[I0, e, theta, polar] = load_itx_LP([foldername,filename]); % Load raw 3d data into I0

% Import the file-tuned parameters from the params table
idx = contains(params.ID, filename(1:-1+find(filename=='_')));
if isempty(idx), idx=1; end
hv = params.hv(idx); 
azi0 = params.azi0(idx); 
tilt0 = params.tilt0(idx); 
pol0 = params.pol0(idx);
FL_shift = params.FLshift(idx);

%% Take mean spectrum edc and set FL as steepest slope
intThPolRng = [];
if isempty(intThPolRng) 
    edc = sum(sum(I0,3),2); % EDC is sum over all polar, theta 
else
    intThIdx = sort( round(interp1( theta, 1:numel(theta), intThPolRng(1,:) )));
    intPolIdx = sort( round(interp1( polar, 1:numel(polar), intThPolRng(2,:) )));
    edc = sum(sum(I0( :, intThIdx(1):intThIdx(2), intPolIdx(1):intPolIdx(2) ),3),2);
end
edcSl = diff(edc); % EDC slope
edcSlx = 0.5*(e(2:end)+e(1:end-1)); % EDC slope x
edcSlxq = edcSlx(1) : 0.001 : edcSlx(end); % Interp EDC slope x
edcSlq = interp1( edcSlx, edcSl, edcSlxq, 'makima'); % Interp EDC 

FL = edcSlxq( find(edcSlq == min(edcSlq))); % Set FL as energy of steepest EDC slope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FL = FL + FL_shift;  %%%%%%%%%% Manually adjust the FL %%%%%%%%
% FL = 111.97;% 97.9806; % Values from SL: 97.9798, 97.9869, 97.9806
disp(['File ',filename,', ','FL = ',num2str(FL),' eV']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure, plot(e, edc,'b'), title({'th,pol int range:';mat2str(intThPolRng);['FL = ',num2str(FL)]});
yyaxis right, plot(edcSlxq, edcSlq,'r'); hold on; 
plot(edcSlx, edcSl, 'm-','LineWidth',1.5);
yLim = ylim; hold on, plot(FL*[1,1], yLim,'k:');
legend({'EDC','Slope','Slope interped','FL'},'Location','southwest');
%% Get constant-energy surface across theta(kx)-vs-polar(ky)

Ecut = 0.00; % Energy to see initial 2d surface
Ehw = 0.01; % Energy half-width to integrate 

thCutVal = 0.43;%0; % Theta val to see a sample spectral cut 
polCutVal = 7;%0; % Polar val to see a sample spectral cut

WF=4.4; V0=10;
geometry = [azi0, tilt0, pol0, hv, WF, V0];
[Esurf0, Kmap0] = constantEsurf(I0, e, theta, polar, FL, 0, Ehw, geometry); % Get initial estimate Fermi Surface 

% Plot constant-energy surface in both degrees, spec cuts at th,pol=0
figure,
subplot(223), imagesc(polar, theta, Esurf0), axis xy
hold on, plot(polCutVal*[1,1],[min(theta),max(theta)],'r--');
hold on, plot([min(polar),max(polar)],thCutVal*[1,1],'r--');
xlabel('polar (deg)'); ylabel('theta (deg)');
title(['Ecut = 0 +/- ',num2str(Ehw),' eV']);

% Calculate index of theta,polar cuts and then plot the spectral cuts 
thCutIdx = round(interp1(theta, 1:numel(theta), thCutVal));
polCutIdx = round(interp1(polar, 1:numel(polar), polCutVal));

subplot(221), imagesc( polar, e-FL, reshape(I0(:,thCutIdx,:), size(I0,1),[]) ); axis xy;
hold on, plot([min(polar),max(polar)],(Ecut-Ehw)*[1,1],'r--');
hold on, plot([min(polar),max(polar)],(Ecut+Ehw)*[1,1],'r--');
title(['th = ',num2str(thCutVal),' deg']);

subplot(224), imagesc( e-FL, theta, I0(:,:,polCutIdx)' ); axis xy; 
hold on, plot((Ecut-Ehw)*[1,1],[min(theta),max(theta)],'r--');
hold on, plot((Ecut+Ehw)*[1,1],[min(theta),max(theta)],'r--');
title(['pol = ',num2str(polCutVal),' deg']);

sgtitle(filename,'Interpreter','none'); colormap turbo;

% Plot 2d surfaces, at several sample energies, in kx-ky
Ecuts = [0, -0.1, -0.5, -1.0];
figure,
for E_i = 1:numel(Ecuts)
    Ecut = Ecuts(E_i);
    [Esurf, Kmap] = constantEsurf(I0, e, theta, polar, FL, Ecut, Ehw, geometry);

    subplot(2,2,E_i);
    pcolor(Kmap(:,:,1), Kmap(:,:,2), Esurf), shading flat;
    hold on, plot([min(min(Kmap(:,:,1))),max(max(Kmap(:,:,1)))],0*[1,1],'w:');
    hold on, plot(0*[1,1],[min(min(Kmap(:,:,2))),max(max(Kmap(:,:,2)))],'w:');
    xlabel('Kx (invA)'); ylabel('Ky (invA)');
    daspect([1,1,1]);
    title(['E = ',num2str(Ecut),' eV'])
end
sgtitle({filename},'Interpreter','none');
colormap turbo

%% Interp 3d I data into an even Kx, Ky grid, Kx-symmetrize, then interp e-axis 

Kmap0x = Kmap0(:,:,1); Kmap0y = Kmap0(:,:,2); % Retrieve Kx,Ky maps at previous E=0 cut, set these as grid limits for all E 

KxLim = min(abs([min(Kmap0x(:)),max(max(Kmap0(:,:,1)))])); % Take 0-centered Kx axis using limits of K map from initial E=0 surface 
Kx = linspace(-KxLim, KxLim, numel(theta)); % theta is in decreasing order, so Kx s/b increasing order
Ky = linspace(max(max(Kmap0(:,:,2))),min(min(Kmap0(:,:,2))), numel(polar)); % polar is increasing order, so Ky s/b decreasing order

[KY,KX] = meshgrid(Ky,Kx); % Rows = Kx, Columns = Ky

a = 1.0*4.321; % Crystal lattice a/x (Angstroms) 
b = 1.0*4.321;  % Crystal lattice b/y

KA = KX/(pi/a); % Plot in BZ-a/x units 
KB = KY/(pi/b); % Plot in BZ-b/y units
 
Iperm = permute(I0, [2,3,1]); % Permute from e,theta,polar to theta,polar,e for per-energy k-calculations
Igrid = zeros( numel(Kx), numel(Ky), numel(e) ); 
for e_i = 1:numel(e)
    EB = (FL-e(e_i));
    Ie = Iperm(:,:,e_i);
    Kmap = theta2kMap( theta, azi0, tilt0, polar+pol0, hv, WF, EB, V0);
    Kmap = permute(Kmap,[2,3,1]); % Permute from Kxyz,theta,pol to theta,pol,Kxyz   
        
    P = [reshape(Kmap(:,:,1),[],1), reshape(Kmap(:,:,2),[],1)]; % Input mat P from scatteredInterpolant documentation 
    F = scatteredInterpolant( P, reshape(Ie,[],1) );
                                
    Ieq = F( KX,KY );
    Igrid(:,:,e_i) = Ieq; 
    
end

% Now interp so e-axis has ~3meV step size 
eq = linspace(e(1), e(end), round(range(e)/0.003) ); 

% Put polar in 3rd dimension to use imresize so 1.E, 2.kx(theta), 3.ky(polar)
Igrid = permute(Igrid, [3,1,2]);  
%%
Iq = imresize( Igrid, [numel(eq), numel(theta)]); % Interpolate to eq
Iq = 0.5 * (Iq + fliplr(Iq)); % Symmetrize across Kx=0 

LV_factor = 0.5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if combine_LH_and_LV == 1
    if contains(filename,'LV'), ILV = Iq; 
    else disp(['Combining with pre-made LV I map...']);    
        I = Iq + LV_factor * ILV; % combine LV + LH 
    end
else, I = Iq;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
FLidx = round(interp1(eq, 1:numel(eq), FL ));
% FS = reshape( I( FLidx, :, : ), numel(Kx), numel(Ky));
FSvol = permute(I,[2,3,1]);
FS = FSvol(:,:,FLidx);

%%%%%%%%%% Symmetrize along ky %%%%%%%%%%%%%%
FS0 = FS; I00 = I;
[BZidxs, FS ] = symmetrize_FSky( FS0, KB );
[~, Iq ] = symmetrize_FSky( I00, KB, [] );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Next step: Evaluate the gap at a given kx,ky coordinate 
eStep = abs(mean(e(2:end)-e(1:end-1))); 
topX = 0.05; % Set energy peak as center-of-mass of top topX % of edc intensity 

KxWidx = round( 0.02 / abs( Kx(2)-Kx(1) )); % Width to integrate Kx 
KyWidx = round( 0.02 / abs( Ky(2)-Ky(1) )); % Width to integrate Ky

gapMap = NaN*ones(numel(Kx),numel(Ky));
for Kx_i = 1:numel(Kx) % increasing Kx may not necessarily be in same direction as I index values
    KxRng = max([1, Kx_i - KxWidx]) : min([numel(Kx), Kx_i + KxWidx]);
    % translate Kx range to I range
    
    for Ky_i = 1:numel(Ky)
        KyRng = max([1, Ky_i - KyWidx]) : min([numel(Ky), Ky_i + KyWidx]);
        spec = nanmean( Iq( :, KxRng, KyRng), 3); 

        [edc, e_centered] = symm_FL_edc(spec, eq, FL );
        
        % Normalize edc for finding peak-com %
        edc = (edc-min(edc))/(max(edc)-min(edc));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        % Gap: find E of peak-com of symmetrized edc 
        % First find peak within 100meV of FL 
        findPkIdx = round(interp1( e_centered, 1:numel(e_centered), [-0.1, 0] )); 
        [pkHt, pkIdx] = max(edc(findPkIdx(1):findPkIdx(2))); 
        pkIdx = pkIdx + findPkIdx-1;
        
        % Find COM using topX% of (normalized) edc intensity 
        topXstart = find( flipud(edc(1:pkIdx)) <= (1-topX)*pkHt, 1, 'first');
        topXstart = pkIdx - topXstart + 1;
        topXend = pkIdx-1 + find( edc(pkIdx:end) <= (1-topX)*pkHt, 1, 'first');        
        topXidx = topXstart : topXend; 
        thePkIdx = round( nansum( topXidx .* edc(topXidx)' )...
                            /nansum(edc(topXidx)) );
                        
        % Full gap is 1x the binding energy of the first peak in FL-symm edc  
        try
            theGap = -1 * abs( e_centered( thePkIdx )); % Using COM pk          
            gapMap( Kx_i, Ky_i ) = theGap;
        end

    end
end

% Plot the gap map and couple Ky spectral cuts 
tlSize = [10,3];%[4,3]; % tiled layout size 
GMpltIdx = 1+[1,2,4,5,7,8,10,11,13,14];%[1,4]; % Position of gap map 
FSpltIdx = 1+[16,17,19,20,22,23,25,26,28,29];% [7,10]; % Position of FS map
subplotIdxs = -2+[3,6;9,12;15,18;21,24;27,30];%[2,5,8,11,3,6,9,12]; % Position of spec cuts 

%
figure; 
% Plot the gap map
subplot(tlSize(1), tlSize(2), GMpltIdx); 
pcolor(KX, KY, gapMap); shading flat; xlabel('Kx'); ylabel('Ky'); title('Gap map'); 
set(gca, 'color',[1,1,1]); colormap((turbo)); daspect([1,1,1]); colorbar();
ax1 = gca; hold(ax1, 'on'); caxis([-0.05,0]);

% Plot the Fermi surface
subplot(tlSize(1), tlSize(2), FSpltIdx); 
pcolor(KX,KY,FS); shading flat; xlabel('Kx'); ylabel('Ky'); title('Fermi surface');
daspect([1,1,1]); colorbar();
ax2 = gca; hold(ax2, 'on');


 % Ky values to take spec cuts 
KyCuts = KyCuts( KyCuts >= min(Ky(:)) & KyCuts <= max(Ky(:)) );
KyCutColors = hsv(numel(KyCuts)); 

MDC_Es = [-0.1 : 0.05 : 0.05]; % Energies to plot MDCs 
EDC_Ks = [-1.0 : 0.2 : 1.0]; % Kx vals to plot EDCs 

MDC_Ehw = round( 0.01 / abs(e(2)-e(1))); % Energy range (idx) to +/- for MDCs
EDC_Khw = round( 0.05 / abs(Kx(2)-Kx(1))); % Kx range (idx) to +/- for EDCs

MDC_Eidxs = round(interp1(e-FL, 1:numel(e), MDC_Es));
EDC_Kidxs = round(interp1(Kx, 1:numel(Kx), EDC_Ks)); 

KyMDCs = zeros(numel(Kx),numel(MDC_Es),numel(KyCuts)); % Array for storing MDCs
KyEDCs = zeros(numel(eq),numel(EDC_Ks),numel(KyCuts)); % Array for storing EDCs

for Ky_i = 1:numel(KyCuts)
    KyCut = KyCuts(Ky_i);
    KyIdx = round( interp1( Ky, 1:numel(Ky), KyCut ));
    
    spec = Iq(:,:,KyIdx);
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
    subplot(tlSize(1),tlSize(2),subplotIdxs(Ky_i,:)); 
    imagesc(Kx, e-FL, spec); axis xy;    ylim([-0.2,0.05]);
    hold on, plot([min(Kx),max(Kx)], 0*[1,1], 'w:');
    title(['Ky = ',num2str(KyCut)]);
    
    % Plot color-coded frame around subplot matching plot on gap map
    boxX = xlim; boxY = ylim; 
    hold on, rectangle('Position',[boxX(1),boxY(1),range(boxX),range(boxY)],'LineWidth',5,'EdgeColor', KyCutColors(Ky_i,:));
    
    if Ky_i ~= numel(KyCuts)
        set(gca,'XTickLabels',{});
    end
    
    % Plot the Ky cut on the gap map 
    plot(ax1, [min(Kx),max(Kx)], KyCut*[1,1], 'Color', KyCutColors(Ky_i,:), 'LineWidth',1.5, 'LineStyle',':')
    plot(ax2, [min(Kx),max(Kx)], KyCut*[1,1], 'Color', KyCutColors(Ky_i,:), 'LineWidth',1.5, 'LineStyle',':')
    
end
sgtitle({filename;['FL = ',num2str(FL),' eV']},'Interpreter','none');

% Plots for Sheng's paper

% Plot the two main panels showing FS, gap map
figure, 
% Plot FS 
subplot(121), 
imagesc(KA(:,1),KB(1,:),FS'); axis xy; 
xlabel('k/(\pi/a)'); ylabel('k/(b/\pi)'); daspect([1,1,1]);
for i = -3:3
    hold on, plot([min(min(KA)),max(max(KA))],i*[1,1],'w-','LineWidth',0.01);
    hold on, plot(i*[1,1],[min(min(KB)),max(max(KB))],'w-','LineWidth',0.01);
end
xlim([-1.8,1.8]); ylim([-3.3,1.3]);
caxis(caxRange(FS,0.05,.995)); colorbar();
    set(gca,'TickDir','out','TickLength',[0.02,.1]);

% Plot the GAP MAP 
subplot(1,2,2), 
imagesc(KA(:,1),KB(1,:),gapMap'); axis xy;
xlabel('k/(\pi/a)'); ylabel('k/(b/\pi)'); daspect([1,1,1]);

for i = -3:3
    hold on, plot([min(min(KA)),max(max(KA))],i*[1,1],'w-','LineWidth',0.01);
    hold on, plot(i*[1,1],[min(min(KB)),max(max(KB))],'w-','LineWidth',0.01);
end
xlim([-1.8,1.8]); ylim([-3.3,1.3]);
caxis([-.035, -0.005]); colorbar();
    set(gca,'TickDir','out','TickLength',[0.02,.1]);


colormap turbo
%%
% Plot panel showing series of deeper-energy FS cuts 
Ecuts = [0,-0.1,-0.5,-1.0];
figure,
for Ei=1:numel(Ecuts)
    [Esurf,~] = constantEsurf(Iq, eq, theta, polar, FL, Ecuts(Ei), Ehw, geometry );
    if Ei==1, Esurf0 = Esurf; end
    
    subplot(1,numel(Ecuts)+1,Ei)
%     pcolor(KA, KB, Esurf); shading flat; 
    imagesc(KA(:,1),KB(1,:), Esurf'); axis xy;
    daspect([1,1,1]);
    xlabel('k/(\pi/a)'); ylabel('k/(b/\pi)'); title(['E = ',num2str(Ecuts(Ei))]);
    for i = -3:3
        hold on, plot([min(min(KA)),max(max(KA))],i*[1,1],'w-','LineWidth',0.1);
        hold on, plot(i*[1,1],[min(min(KB)),max(max(KB))],'w-','LineWidth',0.1);
    end
    caxis(caxRange(Esurf,.2,1))
    xlim([-1.8,1.8]); ylim([-3.3,1.3]);
    set(gca,'TickDir','out','TickLength',[0.02,.1]);
end

EsurfDiff = Esurf0 - Esurf;
% Plot the I difference between E=-1 and E=0 
subplot(1,numel(Ecuts)+1,Ei+1)
imagesc(KA(:,1),KB(1,:), EsurfDiff'); axis xy;
daspect([1,1,1]);
xlabel('k/(\pi/a)'); ylabel('k/(b/\pi)'); 
title(['Diff E=',num2str(Ecuts(1)),' - E=',num2str(Ecuts(Ei))]);% = ',num2str(Ecuts(Ei))]);
for i = -3:3
    hold on, plot([min(min(KA)),max(max(KA))],i*[1,1],'w-','LineWidth',0.1);
    hold on, plot(i*[1,1],[min(min(KB)),max(max(KB))],'w-','LineWidth',0.1);
end
xlim([-1.8,1.8]); ylim([-3.3,1.3]);
colorbar();

colormap turbo

%% Plot slices of the fermi surfaces along z-axis 
[X,Y,Z] = meshgrid(Kx/(pi/a),Ky/(pi/b),eq-FL);
Ip = permute(Iq,[3,2,1]);
Ecutss = {0:.005:.015; 0:-.25:-1};

for Eii= 1:size(Ecutss,1)
    Ecuts = Ecutss{Eii};
    
    figure, 
    slice(X,Y,Z,Ip,[],[],Ecuts); colormap turbo, 
    shading flat, xlabel('Ka'), ylabel('Kb'), zlabel('E (eV)');
    for Ei = 1:numel(Ecuts)
        for i = -2:2
            hold on, plot3([min(min(KA)),max(max(KA))],i*[1,1],Ecuts(Ei)*[1,1],'w-','LineWidth',0.1);
            hold on, plot3(i*[1,1],[min(min(KB)),max(max(KB))],Ecuts(Ei)*[1,1],'w-','LineWidth',0.1);
        end
    end
    caxis([.2,4.2]), pause(.1)
    title(['a=',num2str(a)])
    pause(.1)
    view(45,15)
end

% Symmetrize across kb=0 axis separately for each BZ
%%
[BZidxs, FS_symm2d] = symmetrize_FSky( FS, KB );
BZ0idx = BZidxs{2}; BZ0 = FS_symm2d(:, BZ0idx);
BZ1idx = BZidxs{3}; BZ1 = FS_symm2d(:, BZ1idx);

figure, 
   
% Plot difference between BZ1 BZ0 intensities
BZ10diff = mat2gray(BZ1)-mat2gray(BZ0);
KA10diff = KA(:,BZ0idx); 
KB10diff = KB(:,BZ0idx); 

pcolor(KA10diff, KB10diff, BZ10diff ); axis xy, shading flat
colormap turbo; title(' BZ1 - BZ0 ');
xlim([min(KA10diff(:)),max(KA10diff(:))]);
ylim([min(KB10diff(:)),max(KB10diff(:))]);

% Plot gridlines 
for i = -2:2
    if abs(rem(i,2)) == 1, style = '--'; else, style = ':'; end
    hold on, plot(i*[1,1],[min(KB(:)),max(KB(:))],'w','LineStyle',style);
end
for j = -3:3
    if abs(rem(j,2)) == 1, style = '--'; else, style = ':'; end
        hold on, plot([min(KA(:)),max(KA(:))], j*[1,1],'w','LineStyle',style);
end
daspect([1,1,1]);

% xlim([-1,1]), ylim([-3,1])

end
