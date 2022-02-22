folderName = 'C:\Users\Erica\Desktop\USbTe\2021dec_data\USbTe_1\itx\';
fileExts = {'ust_046_S001.itx'; ...
            'ust_047_S001.itx';...
            'ust_047_S002.itx';...
            'ust_047_S003.itx';...
            'ust_047_S004.itx';...
            'ust_047_S005.itx';...
            'ust_048_S001.itx';...
            'ust_048_S002.itx';...
            'ust_048_S003.itx';...
            'ust_048_S004.itx';...
            'ust_048_S005.itx';...
            'ust_048_S006.itx';...
            'ust_048_S007.itx';...
            'ust_048_S008.itx';...
            'ust_048_S009.itx';...
            'ust_048_S010.itx';...
            'ust_048_S011.itx';...
            'ust_048_S012.itx';...
            'ust_049_S001.itx';...
            'ust_049_S002.itx';...
            'ust_049_S003.itx'};
Zaxis = sort( [134,130:-4:114,106:-7.5:106-7.5*11,17,14,12.6] );
Zaxis_label = 'Temp (K)';
%%
[spec1, eAx, thAx] = load_itx([folderName,fileExts{1}]);
spec3d = zeros(size(spec1,1),size(spec1,2),numel(Zaxis));
spec_i = 1;
for i = numel(fileExts):-1:1
    [spec,eAx,thAx] = load_itx([folderName,fileExts{i}]);
    spec3d(:,:,spec_i) = spec;
    spec_i = spec_i+1;
end
%%
[th, E, Z] = meshgrid( thAx, -eAx, Zaxis);
Zax_int = 12:4:134;
[thq,Eq,Zq] = meshgrid( thAx, -eAx, Zax_int);

spec3d_int = interp3( th,E,Z, spec3d, thq,Eq,Zq, 'spline' );%zeros( size(spec3d,1),size(spec3d,2), numel(Zax_int) );

%%
Ecut = -0.005; 
Ehw = 0.003; % +/- E width
Znorm = 1;
stack = zeros(numel(Zaxis), size(spec3d_int,2));
for i = 1:size(spec3d_int,3)
    Ecut_idx = round(interp1(-eAx, 1:numel(eAx), Ecut));
    Ehw_idx = ceil( Ehw / abs(eAx(2)-eAx(1)) );
    stack_cut = sum( spec3d_int( Ecut_idx-Ehw_idx : Ecut_idx+Ehw_idx, : , i) , 1);   
    if Znorm == 1
        stack_cut = stack_cut / sum(stack_cut) ; 
    elseif Znorm == 2
        stack_cut = stack_cut / sum(sum(spec3d_int(:,:,i)));
    end
    stack(i,:) = stack_cut;
end

figure, 
subplot(223), imagesc(thAx, Zax_int, stack ), axis xy
xlabel('theta (deg)'), 
ylabel(Zaxis_label)
     colormap turbo   

[X,Y] = ginput(1);
Xi = round(interp1( thAx, 1:numel(thAx), X ));
Yi = round(interp1( Zax_int, 1:numel(Zax_int), Y )); 
hold on, plot(thAx(Xi)*[1,1], [Zax_int(1),Zax_int(end)], 'r');
hold on, plot([thAx(1),thAx(end)], Zax_int(Yi)*[1,1], 'r');

Xcut = zeros( size(spec3d_int,3), size(spec3d_int,1));
for i = 1:size(spec3d_int,3)
    Xcut(i,:) = spec3d_int( :, Xi, i ); 
end
Ycut = spec3d_int( :,:,Yi );

subplot(221), imagesc(thAx, -eAx, Ycut), axis xy; 
hold on, plot([thAx(1),thAx(end)], [0,0], 'w:');
title(['T = ',num2str(Zax_int(Yi)),' K']);

subplot(224), imagesc(-eAx, Zax_int, Xcut); axis xy;
hold on, plot([0,0],[Zax_int(1),Zax_int(end)], 'w:');
title(['th = ',num2str(thAx(Xi)),' deg']);

sgtitle(['E = ',num2str(Ecut),' eV']);
%% Sum EDCs to get FL at one or various temperatures (Tidxs) 
edcKrange = round( interp1( thAx, 1:numel(thAx), [-100,100]) );
Tidxs = 1:size(spec3d_int,3);
edcs = cell(1,numel(Tidxs)); 
edcSlopes = cell(1,numel(Tidxs));
FLidxs = zeros(1,numel(Tidxs));

clrs = turbo(numel(Tidxs));
lgd = cell(1,size(clrs,1));
for i = 1:size(spec3d_int,3), lgd{i} = [num2str(Zax_int(i)),'K']; end

for Ti = 1:numel(Tidxs)
    Tidx = Tidxs(Ti);
    edc = sum(spec3d_int(:, max([1,edcKrange(1)]) : min([size(spec3d_int,2),edcKrange(end)]),Tidx),2);
    edcSlope = diff(edc);
    FLidx = 1+find(edcSlope == min(edcSlope));
    edcs{Ti} = edc; 
    edcSlopes{Ti} = edcSlope;
    FLidxs(Ti) = FLidx;
end

figure,
for Ti = 1:numel(Tidxs)
    subplot(221), hold on, plot(eAx, edcs{Ti}, 'Color', clrs(Ti,:))
                hold on, plot(eAx(FLidxs(Ti)), edcs{Ti}(FLidxs(Ti)), 'x','Color',clrs(Ti,:))
                xlabel('E_B (eV)'), ylabel('edc intensity');
    subplot(223), hold on, plot(eAx(2:end), edcSlopes{Ti}, 'Color', clrs(Ti,:))
                hold on, plot(eAx(FLidxs(Ti)), edcSlopes{Ti}(FLidxs(Ti)-1), 'x','Color',clrs(Ti,:))
                xlabel('E_B (eV)'); ylabel('edc slope'); 
    subplot(2,2,[2,4]), hold on, plot(Zax_int(Ti), eAx(FLidxs(Ti)), 'x','Color',clrs(Ti,:),'MarkerSize',10,'LineWidth',2);
        xlabel('T (K)'), ylabel('FL E_B (eV)'); 
end  
legend(lgd)
%% Plot symmetrized-across-FL EDCs versus T 
FLidx = FLidxs(end);
clrs = turbo(size(spec3d_int,3));

for shft = 0%-[.25:.25:1.0]
    edcKrange_deg = [-2.5,-1.5]+shft;%[-6.2,-5.8];%[-.5,.5];%[5.5,6.5];%[-6.5,-5.5];%[-7,-4];%[-1.5,1.5];%[4,7];%
    plt_E_ref = -0.025;

    FLmarg = size(spec3d_int,1) - FLidx;
    symmFLidx = FLidx - FLmarg : FLidx + FLmarg; 
    edcKrange = round( interp1( thAx, 1:numel(thAx), sort(edcKrange_deg) ) );

    edcs = cell(1,size(spec3d_int,3));
    for i = 1:size(spec3d_int,3)
        spec = spec3d_int(:,:,i);

        symmFLspec = spec( FLidx - FLmarg : FLidx + FLmarg, :);

        edc = sum(symmFLspec(:, edcKrange(1):edcKrange(end)), 2);   
        edc = (edc-min(edc))/(max(edc)-min(edc));    
        edc = 0.5*(edc + flipud(edc)); 

        edcs{i} = edc;
    end
    %
    figure,
    for i = 1:size(spec3d_int,3)
        edc = edcs{i} + (i-1)*0.03; 
        hold on, plot(-eAx(FLidx)+eAx( symmFLidx ), edc, 'LineWidth',1.5,'Color',clrs(i,:) );
    end
    yLims = ylim;
    hold on, plot( plt_E_ref*[1,1], [min(yLims),max(yLims)],'k');
    xlim([-.15,.15]);
    title(['k-int [',num2str(thAx(max([1,edcKrange(1)]))),', ',num2str(thAx(min([size(spec3d_int,2),edcKrange(end)]))),']']);
    xlabel(['E_B (eV)']); ylabel('symm. edc int (a.u.)');
end
legend(lgd)