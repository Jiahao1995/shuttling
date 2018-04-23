dotsize = 80;
gapsize = 20:2:40;
tc = zeros(11,41);
barrierHeights = zeros(11,41);
detuningMinimum = zeros(11,41);
% fig = figure;
for ii = 1:11
    for jj = 0:40
        currPotString = [num2str(dotsize) '_' num2str(gapsize(ii)) '_' num2str(jj) '.csv'];
        data = dlmread([sparams.potDir currPotString]);
        [rows,cols] = size(data);
        zdata = data(2:rows,2:cols);
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

        currpot = -interp2(XX,ZZ,zdata,XXq,ZZq);

        % Find which index corresponds to where the 2DEG should be
        [~,index] = min(abs(zz - (-0.5*1E-9)));
        sparams.twoDEGindZ = index;
        
        currTwoDEGPot = currpot(sparams.twoDEGindZ,:);
%         twoDEGPot = squeeze(currpot(:,sparams.twoDEGindZ,:));
        tunBarrierHeight = min(findpeaks(currTwoDEGPot));
        dotMinimum = -max(findpeaks(-currTwoDEGPot));
        barrierHeights(ii,jj+1) = abs(tunBarrierHeight - dotMinimum);
        detuningMinimum(ii,jj+1) = dotMinimum;

        [~,ens] = solve1DSingleElectronSE(sparams,2,xxq*1E-9,currTwoDEGPot*sparams.ee);
        tc(ii,jj+1) = ((ens(2,2) - ens(1,1))/2/sparams.ee);
%         clf;
%         plot(xxq*1E-9,currTwoDEGPot);
%         pause(0.01); 
    end
end
% delete(fig);
barHeights80 = barrierHeights;
detMinimum80 = detuningMinimum;

figure;
volts = linspace(0,2,41);
gap = gapsize;
[VOLTS,GAP] = meshgrid(volts,gap);
surf(VOLTS,GAP,tc);
colormap(jet);
colorbar;
view(2);