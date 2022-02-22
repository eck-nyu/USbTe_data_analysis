
%% With a Kmap and gap map, interpolate into even kx,ky grid and start plotting kx/ky cuts

Kx = Kmap(thCoorsIdx, polCoorsIdx, 1); 
Ky = Kmap(thCoorsIdx, polCoorsIdx, 2);

Kxq = linspace(min(Kx(:)),max(Kx(:)), numel(theta));
Kyq = linspace(min(Ky(:)),max(Ky(:)), numel(polar));

GMq = griddata( reshape(Kx,1,[]), reshape(Ky,1,[]), ...
    reshape(gapMap,1,[]), Kxq, Kyq');

%% There MUST be a way faster way to do this with a single-line function that can take unevenly-spaced 3d-data and interpolate. 
% For now since Kx,Ky are unevenly spaced, just use the griddata and do it for each energy 

Iq = zeros( numel(Kyq), numel(Kxq), numel(e) );
for e_i = 1:numel(e)
    
    Ie = reshape( I(e_i, thCoorsIdx, polCoorsIdx), numel(thCoorsIdx), numel(polCoorsIdx));
    Ieq = griddata( reshape(Kx,1,[]), reshape(Ky,1,[]), reshape(Ie,1,[]), Kxq, Kyq');    
    Iq(:,:,e_i) = Ieq;
end
        
Iqq = permute(Iq, [3,2,1]); 
%% Plot the gap map and couple Ky spectral cuts 
figure; 
subplot(4,3, [1,4,7,10]), 
pcolor(Kxq, Kyq, GMq); shading interp; xlabel('Kx'), ylabel('Ky') ; 
set(gca, 'color',[1,1,1]); colormap turbo; daspect([1,1,1]), 
ax1 = gca; hold(ax1, 'on');

KyCuts = [0, -0.14, -0.30, -0.72, -1.0, -1.22, -1.38]; 
KyCutColors = hsv(numel(KyCuts));
subplotIdxs = [2,5,8,11,3,6,9,12];

for Ky_i = 1:numel(KyCuts)
    KyCut = KyCuts(Ky_i);
    KyIdx = round( interp1( Kyq, 1:numel(Kyq), KyCut ));
    
    % Plot the spectral cuts at each Ky value 
    subplot(4,3,subplotIdxs(Ky_i)); 
    imagesc(Kxq, e-FL, Iqq(:,:,KyIdx)), 
    title(['Ky = ',num2str(KyCut)]);
    axis xy; colormap turbo; ylim([-0.5,0.1]);
    % Plot color-coded frame around subplot matching plot on gap map
    boxX = xlim; boxY = ylim; 
    hold on, rectangle('Position',[boxX(1),boxY(1),range(boxX),range(boxY)],'LineWidth',5,'EdgeColor', KyCutColors(Ky_i,:));
    
    % Plot the Ky cut on the gap map 
    plot(ax1, [min(Kxq),max(Kxq)], KyCut*[1,1], 'Color', KyCutColors(Ky_i,:))
    
end