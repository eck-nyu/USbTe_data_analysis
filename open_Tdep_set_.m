folderName = '/home/data/eck325/USbTe/2021dec_data/USbTe_1/itx/';

% N.E. Photon energy scan [150 : -2 : 92]LH eV
% fileNames = {}; 
% for i = 1:30 
%     fileNumStr = '000'; fileStr = num2str(i); 
%     fileNumStr(end-(numel(fileStr)-1):end)=fileStr; 
%     fileNames{i} = ['ust_036_S',fileNumStr,'.itx'];
% end
% Tvec = [150:-2:92];

% N.E. Photon energy scan [150:-2:92]LV eV
% fileNames = {};
% for i = 1:30
%     fileNumStr = '000'; fileStr = num2str(i);
%     fileNumStr(end-(numel(fileStr)-1):end)=fileStr; 
%     fileNames{i} = ['ust_037_S',fileNumStr,'.itx'];
% end
% Tvec = [150:-2:92];

% The cool-down Tscan @ 112LV, X-G-X 
fileNames = {'ust_046_S001.itx'; ...
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
Tvec = ( [134,130:-4:114,106:-7.5:106-7.5*11,17,14,12.6] );

% % The heating-up Tscan @ 98LH M-X-M
% fileNames = {}; 
% for i = 1:36 
%     fileNumStr = '000'; fileStr = num2str(i); 
%     fileNumStr(end-(numel(fileStr)-1):end)=fileStr; 
%     fileNames{i} = ['ust_039_S',fileNumStr,'.itx'];
% end
% Tvec = linspace(10, 155, numel(fileNames));

[spec3d, E, th, T, FLidx] = open_Tdep_set( folderName, fileNames, Tvec );
if contains(fileNames{1},'046')
    spec3d_112 = spec3d; E_112 = E; th_112 = th; T_112 = T; FLidx_112 = FLidx; 
elseif contains(fileNames{1},'039')
    spec3d_98 = spec3d; E_98 = E; th_98 = th; T_98 = T; FLidx_98 = FLidx; 
end
%%
% For Sheng's paper 
spec3d = spec3d_98;  E = E_98;  th=th_98;  T=T_98;  FLidx = FLidx_98-1; 
    hv = 98; azi0 = .2; tilt0 = 0.2; pol0 = 7;
% spec3d = spec3d_112; E = E_112; th=th_112; T=T_112; FLidx = FLidx_112 +2; 
%     hv = 112; azi0 = .0; tilt0 = .5; pol0 = 0; 
    
% Plot spec
seeTs = sort([10, 50, 80, 120, 135],'descend');
T_specs = cell(1,numel(seeTs)); 
T_kas = zeros(numel(seeTs), numel(th));

EDCcolors = [1,0,0; .5,0,.2; .2,0,.5; .1,0,.7];
see_kas = [0:0.125:1];%[0, 0.3, 1];
ka_w = 12; 
EDCoff = -2; 
EDCxoff = -.1;
EDCyoff = 300;
EDCsig = 3.5;

figure, % Plot spec versus T
for seeT_i = 1:numel(seeTs)
    seeT = seeTs(seeT_i);
    seeTidx = find( abs(T-seeT) == min(abs(T-seeT)) );
    
    V0=10;  
    spec = spec3d(:,:,seeTidx); 
    
    k = theta2kx( th, azi0,tilt0,pol0,hv,WF,0,V0); 
    ka = k / (pi/a); 
    
    % Symmetrize along kx
    [~,k0idx] = min(abs(ka-0));
    flipWidth = min( [k0idx-1, numel(ka)-k0idx] );
    flipIdx = k0idx - flipWidth : k0idx + flipWidth;
    
    spec_symm = spec;
    spec_symm(:,flipIdx) = 0.5* (spec_symm(:,flipIdx) + fliplr(spec_symm(:,flipIdx)));
    
    T_specs{seeT_i} = spec;
    T_kas(seeT_i,:) = ka;
    
%     EDCs = zeros(numel(E), numel(see_kas));
    EDCs = [];
    for see_ka_i = 1:numel(see_kas)
        [~,ka_idx] = min(abs( ka - see_kas(see_ka_i) ));
%         EDC = sum( spec(:, ka_idx-ka_w:ka_idx+ka_w), 2);
% %         EDC = (EDC-min(EDC))/(range(EDC));
%         EDCs(:,see_ka_i) = EDC; 
        
        [EDC,EDC_x] = symm_FL_edc( spec_symm(:, ka_idx-ka_w:ka_idx+ka_w),E,  E(FLidx) );
        EDC = imgaussfilt(EDC,EDCsig);
        EDCs(:,see_ka_i) = EDC;
    end
    EDCs = (EDCs-min(EDCs(:))) / range(EDCs(:));
    EDCs = EDCs - EDCoff*[see_kas];%[0:(numel(see_kas)-1)];
%     EDCs = EDCs -(EDCoff)*see_kas;
%     size(EDC_x')
%     EDC_x = repmat(EDC_x',1,numel(see_kas));% - EDCxoff*see_kas;
%     EDC_x = EDC_x + EDCxoff*[1:numel(see_kas)];
    %     EDC_y = EDCyoff*EDC_x * see_kas'; %repmat( EDCyoff*ones(numel(EDC),1), 1, numel(see_kas));
    EDC_y = EDCyoff*ones(numel(EDC),1) * see_kas;
    
    subplot(numel(seeTs),2,(seeT_i-1)*2+1);
    imagesc(ka, (E-E(FLidx)), spec_symm), axis xy; colormap turbo;        
    
    ylim([-0.145,0.049]);  xlim([-1.1,1.1]);
    caxis(caxRange(spec_symm,.6,1));
    title(['T=',num2str(T(seeTidx))]);
    set(gca,'TickDir','out','TickLength',[0.02,.1]);
    if seeT_i < numel(seeTs), set(gca,'xticklabel',[]); end

    subplot(numel(seeTs),2,(seeT_i-1)*2+2);
    h = plot(EDC_x, EDCs, 'r','LineWidth',1.5);
%     plot3(EDC_x, EDC_y, EDCs)%repmat(see_kas,size(EDC_x,1),1), EDCs,'r');%EDCs, 
	
%     plot(EDC_x, EDCs)
    hold on, plot(0*[1,1],[min(EDCs(:)),max(EDCs(:))],'k:')
    %ylim([min(EDCs(:)),max(EDCs(:))]); 
    xlim([-.145,.145])%[-0.06,.06]);
    set(gca,'TickDir','out','TickLength',[0.02,.1]);
%     xlabel('E'); ylabel('k-cut'); zlabel('EDC int');
%     ylim([0,100])
ylim([min(EDCs(:)),1.1*max(EDCs(:))]);
    if seeT_i < numel(seeTs), set(gca,'xticklabel',[]); end%icks',[]); end
    set(gca,'yticklabel',[]);%icks',[]); 
    set(gca,'ytick',[]);
end




function [spec3d, E, th, T, FLidx] = open_Tdep_set( folderName, fileNames, Tvec )
% folderName, filenames describe location of list of itx spectra to be loaded
% Tdata is vector of z-axis data (usually temperature)

% Load the spectra into 3d array 

% First load 1 sample spec 
[spec, E, th] = load_itx([folderName,fileNames{1}]); % spec 2d array 1.energy (E) 2.theta (th)
if E(end) < E(1), E = -E; end % Flip energy sign if in binding energy

% If energy step is too large, interpolate each spec using imresize to ~3meV energy step 
eStep = abs(E(2)-E(1));
if eStep > 0.005    
    E = E(1) : 0.003 : E(end); % Redefine E to have 3meV step size
    spec = imresize(spec, [numel(E), numel(th)]); % Interped spec
end

spec3d = zeros(size(spec,1),size(spec,2),numel(Tvec)); % 3d array: 1.energy 2.theta 3.temperature
for i = 1:numel(fileNames)
    
    [spec,~,~] = load_itx([folderName,fileNames{i}]); % Assume same E, th for all spec 
    if eStep > 0.005
        spec = imresize( spec, [numel(E), numel(th)] ); 
    end
        
    spec3d(:,:,i) = spec;
end

% Interp T to even steps size, then interp 3d data to interped T values
Tvec_int = linspace(Tvec(1), Tvec(end), round(abs(range(Tvec)/mean(Tvec(2:end)-Tvec(1:end-1)))));

[thGrid, Egrid, Tgrid] = meshgrid( th, E, Tvec);
[thqGrid,EqGrid,TqGrid] = meshgrid( th, E, Tvec_int);


spec3d = interp3( thGrid,Egrid,Tgrid, spec3d, thqGrid,EqGrid,TqGrid, 'makima' );
disp(['Finished loading, interping the data.']);

%%
% specCutCoors = [-0.005, 
if contains(fileNames{1},'046'), E = E+.02; end
Ecuts = [0, -0.01, -0.05, -0.1, -0.2, -0.5]; % Energy value to see iso-E cut 
Ecuts = [0];
Ehw = 0.003; % +/- E width

for Ecut_i = 1:numel(Ecuts)
    Ecut = Ecuts(Ecut_i);

    EcutIdx = round(interp1(E, 1:numel(E), Ecut));
    EhwIdx = ceil( Ehw / abs(E(2)-E(1)) );

    normInt = 1; % Normalize over T by 1: set sum MDC to 1, 2: set sum of entire spec to 1
    EcutSpec = zeros(numel(Tvec), size(spec3d,2));
    for i = 1:size(spec3d,3)
        spec = spec3d(:,:,i);

        Ecut_th = sum( spec( EcutIdx-EhwIdx : EcutIdx+EhwIdx, : ) , 1);   
        if normInt == 1
            Ecut_th = Ecut_th / sum(Ecut_th) ; 
        elseif normInt == 2
            Ecut_th = Ecut_th / sum(sum(spec));
        end
        EcutSpec(i,:) = Ecut_th;
    end
    if Ecut==0, EcutSpec0 = EcutSpec; end
        
    figure, 
    subplot(223), imagesc(th, Tvec_int, EcutSpec ), axis xy
    xlabel('theta (deg)'), ylabel('T (K)')
    title(['E = ',num2str(Ecut),' eV']);

    % Get click input for horizontal, vertical spec cuts
    [thCut,Tcut] = ginput(1);
    thIdx = round(interp1( th, 1:numel(th), thCut ));
    Tidx = round(interp1( Tvec_int, 1:numel(Tvec_int), Tcut )); 

    hold on, plot(th(thIdx)*[1,1], [Tvec_int(1),Tvec_int(end)], 'r');
    hold on, plot([th(1),th(end)], Tvec_int(Tidx)*[1,1], 'r');

    thSpec = zeros( size(spec3d,3), size(spec3d,1), numel(thCut));
    for i = 1:size(spec3d,3)
        thSpec(i,:) = spec3d( :, thIdx, i ); 
    end
    Ycut = spec3d( :,:,Tidx );

    % Plot the E-theta cut
    subplot(221), imagesc(th, E, Ycut), axis xy; 
    hold on, plot([th(1),th(end)], [0,0], 'w:');
    title(['T = ',num2str(Tvec_int(Tidx)),' K']);
    ylim([-0.2, 0.05]); ylabel('E (eV)'); xlabel('theta (deg)');
    

    % Plot the theta-T cut 
    subplot(224), imagesc(E, Tvec_int, thSpec); axis xy;
    hold on, plot([0,0],[Tvec_int(1),Tvec_int(end)], 'w:');
    xlim([-0.2, 0.05]); xlabel('E (eV)'); ylabel('T (K)');
    title(['th = ',num2str(th(thIdx)),' deg']);
    
    colormap turbo
end
pause(.1);

%% Sum EDCs to get FL at one or various temperatures (Tidxs), over 

figure,  
subplot(223), imagesc(th, Tvec_int, EcutSpec0 ), axis xy, colormap turbo;
xlabel('theta (deg)'), ylabel('T (K)')
title(['E = 0 eV']);

sgtitle({['Select theta range to sum over for edc to set FL,'];['T range to average over FL to set gap']});

% Get click input for th-range to see FL(T) and gap(T)
[thCut,Tcut] = ginput(2);
hold on, rectangle('Position',[min(thCut),min(Tcut),range(thCut),range(Tcut)],'EdgeColor','w');

edcThrange = round( interp1( th, 1:numel(th), thCut') );
Tidxs = 1:size(spec3d,3);
edcs = cell(1,numel(Tidxs)); 
edcSlopes = cell(1,numel(Tidxs));
FLidxs = zeros(1,numel(Tidxs));

clrs = turbo(numel(Tidxs));
lgd = cell(1,size(clrs,1));
for i = 1:size(spec3d,3), lgd{i} = [num2str(Tvec_int(i)),'K']; end

for Ti = 1:numel(Tidxs)
    Tidx = Tidxs(Ti);
    edc = sum(spec3d(:, max([1,edcThrange(1)]) : min([size(spec3d,2),edcThrange(end)]),Tidx),2);
    edcSlope = diff(edc);
    FLidx = 1+find(edcSlope == min(edcSlope));
        FLidx = FLidx(1); 
    edcs{Ti} = edc; 
    edcSlopes{Ti} = edcSlope;
    FLidxs(Ti) = FLidx;
end

% Plot the FL versus T summary 
figure,
for Ti = 1:numel(Tidxs)
    subplot(221), hold on, plot(E, edcs{Ti}, 'Color', clrs(Ti,:))
                hold on, plot(E(FLidxs(Ti)), edcs{Ti}(FLidxs(Ti)), 'x','Color',clrs(Ti,:))
                xlabel('E_B (eV)'), ylabel('edc intensity');
                title('EDC over T');
    subplot(223), hold on, plot(E(2:end), edcSlopes{Ti}, 'Color', clrs(Ti,:))
                hold on, plot(E(FLidxs(Ti)), edcSlopes{Ti}(FLidxs(Ti)-1), 'x','Color',clrs(Ti,:))
                xlabel('E_B (eV)'); ylabel('edc slope'); 
                title('EDC slope over T');
    subplot(2,2,[2,4]), hold on, plot(Tvec_int(Ti), E(FLidxs(Ti)), 'x','Color',clrs(Ti,:),'MarkerSize',10,'LineWidth',4);
        xlabel('T (K)'), ylabel('FL E_B (eV)');  
        title('FL versus T');
end  
%% Plot gap versus T 
FLidx = round(mean(FLidxs)); % FLidxs(end);
clrs = turbo(size(spec3d,3));

for shft = 0
    edcKrange_deg = thCut' + shft; %[-2.5,-1.5]+shft;
    plt_E_ref = -0.025;

    FLmarg = size(spec3d,1) - FLidx;
    symmFLidx = FLidx - FLmarg : FLidx + FLmarg; 
    edcThrange = round( interp1( th, 1:numel(th), sort(edcKrange_deg) ) );

    edcs = cell(1,size(spec3d,3));
    for i = 1:size(spec3d,3)
        spec = spec3d(:,:,i);

        symmFLspec = spec( FLidx - FLmarg : FLidx + FLmarg, :);

        edc = sum(symmFLspec(:, edcThrange(1):edcThrange(end)), 2);   
        edc = (edc-min(edc))/(max(edc)-min(edc));    
        edc = 0.5*(edc + flipud(edc));

        edcs{i} = edc;
    end
    %
    figure,
    for i = 1:size(spec3d,3)
        edc = edcs{i} + (i-1)*0.03; 
        hold on, plot(E(FLidx)-E( symmFLidx ), edc, 'LineWidth',1.5,'Color',clrs(i,:) );
    end
    yLims = ylim;
    hold on, plot( plt_E_ref*[1,1], [min(yLims),max(yLims)],'k');
    xlim([-.15,.15]);
    title(['k-int [',num2str(th(max([1,edcThrange(1)]))),', ',num2str(th(min([size(spec3d,2),edcThrange(end)]))),']']);
    xlabel(['E_B (eV)']); ylabel('symm. edc int (a.u.)');
end
legend(lgd)

T = Tvec_int;
end