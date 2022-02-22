function [edc, e_centered] = symm_FL_edc( spec, e, FL )
% Sums spectra over columns for edc, symmetrizes across FL (if only
% considering sub-range of theta/kx, crop spec before input)

% Order: 1. Get EDC, 2. Normalize EDC, 3. FL-symmetrize last

FLidx = round(interp1( e, 1:numel(e), FL ));
FLmargin = size(spec,1) - FLidx;

% Crop lower energies so there's equal pixels above/below FL 
symmIdx = FLidx - FLmargin : FLidx + FLmargin;

edc = sum( spec(symmIdx, :), 2);
%edc = (edc-min(edc))/(max(edc)-min(edc));
edc = 0.5 * (edc + flipud(edc));

% e_cropped = e( symmIdx );
e_centered = e( symmIdx ) - FL;
end