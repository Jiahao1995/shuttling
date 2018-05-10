shuttleParameterFile;

% dotsize = 40:2:80;
% gapsize = 20:2:40;
% vl = 0.1:0.1:2.0;

% dotsize = 20:10:80;
dotsize = [40,80];
gapsize = 5:5:50;
vl = 0.5:0.1:2.5;

dir = 'simulatedPotentials/Two Gates with the Same Voltages/';

pp = 1;
% figure;
for ii = 1:length(dotsize)
    fprintf(1, '\nCurrently on dot size %d\n', dotsize(ii));
    tc = zeros(length(dotsize), length(vl), length(gapsize));
    difference = zeros(length(vl), length(gapsize));
    threshold = 0.001;
    for jj = 1:length(gapsize)
        for kk = 1:length(vl)
            currPotString = ['DOT_SIZE_' num2str(dotsize(ii)) '_GAP_SIZE_' num2str(gapsize(jj)) '_VOLTAGE_' sprintf('%.1f', vl(kk)) '.csv'];
            data = dlmread([dir currPotString]);
            [rows, cols] = size(data);
            zdata = data(2:rows, 2:cols);
            xdata = data(1, 2:cols);
            ydata = data(2:rows, 1);
            
            desiredGridX = 2^(nextpow2(length(xdata)));
            desiredGridZ = 3*length(ydata);
            
            xxq = linspace(min(xdata), max(xdata), desiredGridX);
            zzq = linspace(min(ydata), max(ydata), desiredGridZ);
            
            [XX, ZZ] = meshgrid(xdata, ydata);
            [XXq, ZZq] = meshgrid(xxq, zzq);
            
            currpot = -interp2(XX, ZZ, zdata, XXq, ZZq);
            [~, index] = min(abs(zzq - (-0.5)));
            twoDEGIndZ = index;
            
            currTwoDEGPot = currpot(twoDEGIndZ, :);
            

            max_potential = findpeaks(currTwoDEGPot);
            min_potential = -max(findpeaks(-currTwoDEGPot));
            if isempty(max_potential)
                max_potential = min_potential;
            end
            difference(kk, jj) = max_potential - min_potential;
            if difference(kk,jj) < threshold
                difference(kk,jj) = NaN;
            end
   
            [~, ens] = solve1DSingleElectronSE(sparams, 2, xxq*1E-9, currTwoDEGPot*sparams.ee);
            deltaE = ens(2,2) - ens(1,1);            
            tc(ii, kk, jj) = 0.5 * abs(deltaE);
            if isnan(difference(kk,jj))
                tc(ii,kk,jj) = NaN;
            end
        end 
    end
%     figure;
%     [VL, GAP] = meshgrid(vl, gapsize);
% %     size(VL)
% %     size(GAP)
% %     size(tc)
%     s = surf(VL, GAP, tc(:,:)'/sparams.ee*1E3);
%     view(2);
%     ylabel('Gap Size [nm]');
%     xlabel('Potential [V]');
%     set(s, 'edgecolor', 'none');
%     xlim([min(vl), max(vl)]);
%     title(sprintf('Tunnel Coupling t_c with Dot Size %d', dotsize(ii)));
%     colormap(jet);
%     colorbar;
%     drawnow;
    pp = pp + 1;

    if mod(ii-1,6) == 0
        fig = figure;
        pos = get(fig,'position');
        set(fig,'position',[pos(1:2)/4 pos(3)*2.0 pos(4)*1.25]);
    end
    subplot(2,3,mod(ii-1,6)+1);
    [VOLTS,GAP] = meshgrid(vl,gapsize);
    
    %if ii == 1
    %    lims = [min(min(difference)),max(max(difference))];
    %end
    %s = surf(VOLTS,GAP,difference');
%     
%     if ii == 1
%         lims = [min(min(min(tc))),max(max(max(tc)))];
%     end
    s = surf(VOLTS,GAP,squeeze(tc(ii,:,:))');
    
    title(sprintf('Dot Size = %d [nm]',dotsize(ii)));
    xlabel('Voltage [V]');
    ylabel('Gap Size [nm]');
    set(s,'edgecolor','none');
    xlim([min(min(VOLTS)),max(max(VOLTS))]);
    ylim([min(min(GAP)),max(max(GAP))]);
    caxis(lims);
    colormap(jet);
    colorbar;
    view(2);
    drawnow;
end


%%
% Plot individual functions
dsize = 80;
gsize = 15;
v = 0.5;
currPotString = ['DOT_SIZE_' num2str(dsize) '_GAP_SIZE_' num2str(gsize) '_VOLTAGE_' sprintf('%.1f', v) '.csv'];
data = dlmread([dir currPotString]);
[rows, cols] = size(data);
zdata = data(2:rows, 2:cols);
xdata = data(1, 2:cols);
ydata = data(2:rows, 1);

desiredGridX = 2^(nextpow2(length(xdata)));
desiredGridZ = 3*length(ydata);

xxq = linspace(min(xdata), max(xdata), desiredGridX);
zzq = linspace(min(ydata), max(ydata), desiredGridZ);

[XX, ZZ] = meshgrid(xdata, ydata);
[XXq, ZZq] = meshgrid(xxq, zzq);

currpot = -interp2(XX, ZZ, zdata, XXq, ZZq);
[~, index] = min(abs(zzq - (-0.5)));
twoDEGIndZ = index;

currTwoDEGPot = currpot(twoDEGIndZ, :);

[wfs, ens] = solve1DSingleElectronSE(sparams, 2, xxq*1E-9, currTwoDEGPot*sparams.ee);

figure;
hold on;
plot(xxq,currTwoDEGPot);
plot(xxq,wfs.^2/500000000 + min(currTwoDEGPot));


