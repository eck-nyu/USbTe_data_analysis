function [ spec3d, Eaxis, angleAxis, polarAxis ] = load_itx_LP( fileName )
% Taking Yishuai's load_itx function, tweaking to take Jonathan's LP scans
% (Live Polar) which are sets of spectral cuts (energy-vs-theta) taken over
% range of polar angles. Output is 3d image array, with the energy, theta, 
% and polar axes 
    
    fdata=importdata(fileName);
    
    spec=fdata.data;
    
    % now load the energy and angle axis
    
    f1=fopen(fileName);
    
    for i=1:2
        fline=fgetl(f1);
        
        % now import data size
        t_size=cell2mat(textscan(fline,'WAVES/D/N=(%f, %f, %f)'));
    end
    
    while ~strcmp(fline,'END')
        fline=fgetl(f1);
    end
    % now read in the next line contain the energy and angle axies
    fline=fgetl(f1);
    
    % split the line
    line_string=split(fline(2:end),';');
    xinfo=cell2mat(textscan(line_string{1},' SetScale/P x %f, %f'));
    yinfo=cell2mat(textscan(line_string{2},' SetScale/P y %f, %f'));
    
    %%%% 2022-01-11 Added these two lines for opening LivePolar scans from JD %%%% 
    zinfo = cell2mat(textscan(line_string{3},' SetScale/P z %f, %f'));
    polarAxis = zinfo(1)+[0:t_size(3)-1]*zinfo(2);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
    Eaxis=xinfo(1)+[0:t_size(1)-1]*xinfo(2);
    angleAxis=yinfo(1)+[0:t_size(2)-1]*yinfo(2);

    spec=spec(:,~all(isnan(spec)));
    % close file
    fclose(f1);
 
    %%% 2022-01-13 Added lines to reshape, permute data into conventional
    %%% 3d array %%%%%%%%%%%%%    
    % Reshape the data to 3d array
    spec3d = zeros( numel(Eaxis), numel(polarAxis), numel(angleAxis));
    for i = 1:numel(angleAxis)
        theSpec = reshape(spec(:,i), numel(Eaxis), numel(polarAxis));
        spec3d(:,:,i) = theSpec;
    end
    % Permute axes so 1=energy, 2=theta, 3=polar 
    spec3d = permute(spec3d,[1,3,2]);     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

