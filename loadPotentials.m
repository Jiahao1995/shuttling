function [ xx, zz, pots ] = loadPotentials( sparams )
%LOADPOTENTIALS Function to either generate the potential profiles
%automatically or load them from an external file
    
    nn = 1;
    for ii = 1:4
        for jj = 0:59
            currPotString = ['P' num2str(ii) '_' num2str(jj) '.csv'];
            
            % Load the figure
            data = dlmread([sparams.potDir currPotString]);
            [rows,cols] = size(data);
            zdata = data(2:rows,2:cols);
            
            % Then we are on the first potential
            if nn == 1
                xdata = data(1,2:cols);
                ydata = data(2:rows,1);
                
                % The data may not be uniform in sampling, so we need to fix that for
                % the fourier transforms in the main code for speed up.
                % Find next highest power of 2 to the length of xx
                desiredGridX = 2^(nextpow2(length(xdata)));
                desiredGridZ = 3*length(ydata);
                
                % Make linearly spaced grid of points to interpolate the
                % potential at
                xxq = linspace(min(xdata),max(xdata),desiredGridX);
                zzq = linspace(min(ydata),max(ydata),desiredGridZ);
                
                [XX,ZZ] = meshgrid(xdata,ydata);
                [XXq,ZZq] = meshgrid(xxq,zzq);
                
                pots = zeros(1,desiredGridZ,desiredGridX);
            end
            pots(nn,:,:) = interp2(XX,ZZ,zdata,XXq,ZZq);
            nn = nn + 1;
        end
    end
    pots = -pots*sparams.ee; % Convert to J
    xx = xxq*1E-9; % Convert to m
    zz = zzq*1E-9; 
end

