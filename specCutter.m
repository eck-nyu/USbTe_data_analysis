function specCut = specCutter(Idata, e, theta, polar, Kmap, axisPts, axisLabels, axisNumPts)

axisVec = [ linspace(axisPts(1,1), axisPts(2,1), axisNumPts)',...
            linspace(axisPts(1,2), axisPts(2,2), axisNumPts)'];

        
specCut = zeros(numel(e), axisNumPts);

for ax_i = 1:axisNumpts
    x = axisVec(ax_i,1);
    y = axisVec(ax_i,2);
    
    ptDist = sqrt( (Kmap(:,:,1)-x).^2 + (Kmap(:,:,2)-y).^2 );
    [rowIdx,colIdx] = find(ptDist == nanmin(ptDist(:)),1,'first'); 

%     specCut(:, ax_i) = Idata(:, 


end