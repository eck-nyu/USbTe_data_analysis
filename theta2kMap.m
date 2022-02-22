
function kMap = theta2kMap( thetaVec, azi0, tilt0, polarVec, hv, WF, EB, innerPot )
    % Enter all angles in degrees, all energies in eV
    % Edited 2022/01/27 to incorporate vector of pol0 values 
    
    ec = 1.6e-19; me=9.109e-31; hbar=1.0546e-34; 
    
    kMap = zeros(3, numel(thetaVec), numel(polarVec)); % Empty array for whole theta-polar space  
    for pol0_i = 1:numel(polarVec)
        pol0 = polarVec(pol0_i); 
        
        % Create the rotation matrix for set of angles 
        eul = [azi0, tilt0, pol0]*pi/180; 
        rotm = eul2rotm(eul,'ZYX');
        
        kVec = zeros(3, numel(thetaVec)); % Starting vecs for kx, ky, kz for each polar coordinate
        
        kVec(1,:) = sin(thetaVec'*pi/180); % Init kx is sin of thetax
        kVec(2,:) = zeros(numel(thetaVec),1); % Init ky is 0
        kVec(3,:) = cos(thetaVec'*pi/180); % Init kz is cos of thetax
        
        kVec = (rotm*kVec); % Apply the total rotation matrix to kMap. 

        kVec(1:2,:)=kVec(1:2,:)*(sqrt(((hv-WF-EB)*ec)*2*me)/hbar)*1e-10;
        kVec(3,:)=(sqrt((((hv-WF-EB)*kVec(3,:).*kVec(3,:)+innerPot)*ec)*2*me)/hbar)*1e-10;
        
        kMap(:,:,pol0_i) = kVec;
    end
end

