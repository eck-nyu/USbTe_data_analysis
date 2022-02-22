function [BZidxs, FS_symm] = symmetrize_FSky( FS, KB, Enorm )
% FS 2 or 3d matrix with 1.ky 2.kx (3.E)

% Unfortunately a 2d FS dimensions are 1.kx, 2.ky
% and a 3d I dimensions are 1.E, 2.kx, 3.ky

% If given a 3d I, permute so E is 3rd dimension
if size(FS,3) > 1 
    FS = permute( FS, [2,3,1] ); % 1.kx, 2.ky, 3.E
end

% Find integer kb(ky) idx for each BZ 
% [~,BZ0idx_start] = min( abs( KB(1,:) - 1 ));
% [~,BZ0idx_end]   = min(abs( KB(1,:) + 1 ));
% [~,BZ1idx_end]   = min(abs( KB(1,:) + 3 ));
% 
% % Set idx ranges for each BZ 
% BZm1idx = 1:BZ0idx_start-1; 
% BZ0idx = BZ0idx_start : BZ0idx_end; 
% BZ1idx = BZ0idx_end + 1 : BZ1idx_end + 1;
% BZ2idx = BZ1idx_end + 2 : size(FS,2);

% Another attempt to define BZs 
[~,BZ0idx_center] = min( abs( KB(1,:) ));
[~,BZ1idx_center] = min( abs( KB(1,:) + 2));
BZhw = round( 0.5 * abs( BZ1idx_center - BZ0idx_center ) );

BZm1idx = 1:BZ0idx_center-BZhw;
BZ0idx = BZ0idx_center + (-BZhw:BZhw);
BZ1idx = BZ1idx_center + (-BZhw:BZhw);
BZ2idx = BZ1idx_center+BZhw : size(FS,2);


% Get each BZ intensity from FS 
BZm1 = FS(:, BZm1idx, :);
BZ0 = FS( :, BZ0idx, : ); 
BZ1 = FS( :, BZ1idx, : );
BZ2 = FS( :, BZ2idx, : ); 

BZm1 = 0.5*( BZm1 + fliplr(BZ1(1:size(BZm1,2))) );
BZ1 = 0.5*( BZ1 + fliplr(BZ1) );
BZ2 = 0.5*( BZ2 + fliplr(BZ0(end-size(BZ2,2)+1:end) ) );
BZ0 = 0.5*( BZ0 + fliplr(BZ0) );

% if renormI == 1
if size(FS,3) > 1   
    BZm1 = BZm1 / mean(mean(BZm1(:, end,:))) * mean(mean(BZ0(:,1,:)));  %(BZ1m - min(BZ1m(:))) / range(BZ1m(:));
    BZ1 = BZ1 / mean(mean(BZ1(:,1,:))) * mean(mean(BZ0(:,end,:))); %mean(BZ(BZ1 - min(BZ1(:))) / range(BZ1(:)); %BZ1 / sum(BZ1(:)) * size(BZ1,1);
    BZ2 = BZ2 / mean(mean(BZ2(:,1,:))) * mean(mean(BZ1(:,end,:)));%(BZ2 - min(BZ2(:))) / range(BZ2(:)); %BZ2 / sum(BZ2(:));% * size(BZ2,1);

else 
    % Renormalize intensities per BZ at BZ boundaries 
    BZm1 = BZm1 / mean(BZm1(:, end)) * mean(BZ0(:,1));  
    BZ1 = BZ1 / mean(BZ1(:,1)) * mean(BZ0(:,end)); 
    BZ2 = BZ2 / mean(BZ2(:,1)) * mean(BZ1(:,end));

end

% Put the BZs together into symm FS 
FS_symm = NaN*ones(size(FS));
% FS_symm(:, 1:size(BZ1m,2), :) = BZ1m;
% FS_symm(:, BZ0idx_start-1 + (1:size(BZ0,2)), : ) = BZ0;
% FS_symm(:, BZ0idx_start-1+size(BZ0,2) + (1:size(BZ1,2)), : ) = BZ1; 
% FS_symm(:, end-size(BZ2,2)+1:end, : ) = BZ2;

FS_symm(:, BZm1idx, :) = BZm1;
FS_symm(:, BZ0idx, :) = BZ0;
FS_symm(:, BZ1idx, :) = BZ1;
FS_symm(:, BZ2idx, :) = BZ2;

BZidxs = {BZm1idx; BZ0idx; BZ1idx; BZ2idx};

if size(FS,3) > 1
    FS_symm = permute(FS_symm, [3,1,2]); % Back to 1.E, 2.kx, 3.ky
end

end




% function [BZidxs, FS_symm] = symmetrize_FSky( FS, KB )
% 
% % Find integer kb(ky) idx for each BZ 
% 
% [~,BZ0idx_start] = min( abs( KB(1,:) - 1 ));
% [~,BZ0idx_end]   = min(abs( KB(1,:) + 1 ));
% [~,BZ1idx_end]   = min(abs( KB(1,:) + 3 ));
% 
% % Set idx ranges for each BZ 
% BZ1midx = 1:BZ0idx_start-1; 
% BZ0idx = BZ0idx_start : BZ0idx_end; 
% BZ1idx = BZ0idx_end + 1 : BZ1idx_end + 1;
% BZ2idx = BZ1idx_end + 2 : size(FS,2);
% 
% % Get each BZ intensity from FS 
% BZ1m = FS(:, BZ1midx );
% BZ0 = FS( :, BZ0idx ); 
% BZ1 = FS( :, BZ1idx );
% BZ2 = FS( :, BZ2idx ); 
% 
% size(BZ0),size(BZ1)
% 
% BZ1m = 0.5*( BZ1m + fliplr(BZ1(1:size(BZ1m,2))) );
% BZ1 = 0.5*( BZ1 + fliplr(BZ1) );
% BZ2 = 0.5*( BZ2 + fliplr(BZ0(end-size(BZ2,2)+1:end) ) );
% BZ0 = 0.5*( BZ0 + fliplr(BZ0) );
% 
% BZ1m = BZ1m / mean(BZ1m(:, end)) * mean(BZ0(:,1));  %(BZ1m - min(BZ1m(:))) / range(BZ1m(:));
% BZ1 = BZ1 / mean(BZ1(:,1)) * mean(BZ0(:,end)); %mean(BZ(BZ1 - min(BZ1(:))) / range(BZ1(:)); %BZ1 / sum(BZ1(:)) * size(BZ1,1);
% BZ2 = BZ2 / mean(BZ2(:,1)) * mean(BZ1(:,end));%(BZ2 - min(BZ2(:))) / range(BZ2(:)); %BZ2 / sum(BZ2(:));% * size(BZ2,1);
% 
% FS_symm = NaN*ones(size(FS));
% FS_symm(:, 1:size(BZ1m,2) ) = BZ1m;
% FS_symm(:, BZ0idx_start-1 + (1:size(BZ0,2))) = BZ0;
% FS_symm(:, BZ0idx_start-1+size(BZ0,2) + (1:size(BZ1,2)) ) = BZ1; 
% FS_symm(:, end-size(BZ2,2)+1:end) = BZ2;
% 
% BZidxs = {BZ1midx; BZ0idx; BZ1idx; BZ2idx};
% end
