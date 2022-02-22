% Try to get a spectral cut along specified axis
% axisPts 1st row: starting kx,ky coordinate; 2nd row, ending 
axisPts = [0.7,-0.7; 0.0, 0.0]; 
axisLabels = {'M','\Gamma'};

axisPts = [0.0, 0.0; 0.7,0]; 
axisLabels = {'\Gamma','X'};

axisPts = [0.7, 0; 0.7, -0.7];
axisLabels = {'X','M'};

axisNumPts = 80; % Number of pts to sample along axis 
axisVec = [ linspace(axisPts(1,1), axisPts(2,1), axisNumPts)',...
            linspace(axisPts(1,2), axisPts(2,2), axisNumPts)'];

Iperm = permute(I0,[2,3,1]); 
axisCut = zeros( numel(e), axisNumPts );
for e_i = 1:numel(e)
    Ie = Iperm(:,:,e_i);
    Kmap = theta2kMap( theta, azi0, tilt0, polar+pol0, hv, WF, EB, V0);
    Kmap = permute(Kmap,[2,3,1]); % Permute from Kxyz,theta,pol to theta,pol,Kxyz
    P = [reshape(Kmap(:,:,1),[],1), reshape(Kmap(:,:,2),[],1)]; % Input mat P from scatteredInterpolant documentation 
    F = scatteredInterpolant( P, reshape(Ie,[],1) );
       
    for ax_i = 1:axisNumPts
        
        Ival = F( axisVec(ax_i,1), axisVec(ax_i,2) );
        axisCut(e_i, ax_i) = Ival;
    end
end
axisCutq = imresize(axisCut, [numel(eq), axisNumPts]);
%%
figure, imagesc([1:axisNumPts],eq-FL,axisCutq), axis xy
hold on, plot([1,axisNumPts], 0*[1,1], 'w:'); 
xticks([1,axisNumPts]); xticklabels(axisLabels)
ylim([-0.2,0.05]);
colormap turbo