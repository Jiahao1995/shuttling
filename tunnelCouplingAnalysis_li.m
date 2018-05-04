shuttleParameterFile;
% dotsize = 35:5:55;
% gapsize = 20:5:40;
% dotsize = 35;
% dotsize = 20:10:30;
% dotsize = 20:10:30;
% gapsize = 16:8:32;
dotsize = 30;
gapsize = 24;
% gapsize = 20:5:40;
% gapsize = 20:5:40;
% vl = [0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9];
% vl = [1.0,1.5,2.0];
% vl = [0.2,0.4:0.2:1.9];
% vr = 0:1:10;
% vl = 0.6:0.1:1.9;
vl = 0.6:0.1:1.4;
% vl = 0.6;
% vr = 0:10;
% vr = 0:10;
diff = [-0.2,-0.192,-0.184,-0.176,-0.168,-0.16,-0.152,00.144,-0.136,-0.128,-0.12,-0.112,-0.104,-0.096,-0.088,-0.08,-0.072,-0.064,-0.056,-0.048,-0.04,-0.032,-0.024,-0.016,-0.008,0,0.008,0.016,0.024,0.032,0.04,0.048,0.056,0.064,0.072,0.08,0.088,0.096,0.104,0.112,0.12,0.128,0.136,0.144,0.152,0.16,0.168,0.176,0.184,0.192,0.2];
% diff = linspace(-0.003,0.003,41);
dir = 'simulatedPotentials/0427/';

pp = 1;
% figure;
for ii = 1:length(dotsize)
    for jj = 1:length(gapsize)
        fprintf(1,'\nCurrently on dot size %d and gap size %d\n',dotsize(ii), gapsize(jj));
        tc = zeros(length(vl),length(diff));
        tcOther = zeros(length(vl),length(diff));
        for kk = 1:length(vl)
            for ll = 1:length(diff)
%                 currPotString = [num2str(dotsize(ii)) '_' num2str(gapsize(jj)) '_' num2str(vl(kk)) '_' num2str(vr(ll)) '.csv'];
                % currPotString = [num2str(dotsize(ii)) '_' num2str(gapsize(jj)) '_' sprintf('%.1f',vl(kk)) '_' num2str(vr(ll)) '.csv'];
                currPotString = ['DOT_SIZE_' num2str(dotsize(ii)) '_GAP_SIZE_' num2str(gapsize(jj)) '_VL_' sprintf('%.1f',vl(kk)) '_DIFF_' num2str(diff(ll)) '.csv'];
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
                
                currTwoDEGPot = currpot(twoDEGIndZ,:); 
%                 findpeaks(currTwoDEGPot);
                [~, ens] = solve1DSingleElectronSE(sparams, 2, xxq*1E-9, currTwoDEGPot*sparams.ee);
                deltaE_ = ens(2,2) - ens(1,1);
%                 
                peaks = findpeaks(-currTwoDEGPot);
%                 figure;
%                 plot(xxq,currTwoDEGPot);
% %                 delay(0.3);
%                 pause(0.15);
%                 clf;
                min_potential = min(peaks);
                max_potential = max(peaks);
                deltaE = (max_potential - min_potential) * sparams.ee;
                
                if deltaE_ < deltaE
                    tcOther(kk, ll) = NaN;
                else
                    tcOther(kk, ll) = sqrt(deltaE_^2 - deltaE^2)/2;
                end               

                tc(kk,ll) = calculateTunnelCoupling(sparams, xxq*1E-9, currTwoDEGPot*sparams.ee);
            end
        end
%         if mod(pp-1,6) == 0
%             fig = figure;
%             pos = get(fig,'position');
%             set(fig,'position',[pos(1:2)/4 pos(3)*2.0 pos(4)*1.25]);
%         end
%         subplot(2,3,mod(pp-1,6)+1);
%         subplot_index = (ii-1) * length(dotsize) + jj;
%         subplot(length(dotsize), length(gapsize), subplot_index)
        figure;
        [DIFF, VL] = meshgrid(diff, vl);
        s = surf(DIFF, VL, tc(:,:)/sparams.ee*1E3);
        ylabel('V_L [V]');
        xlabel('V_L-V_R [V]');
        set(s,'edgecolor','none');
        view(2);
        xlim([min(diff),max(diff)]);
        title(sprintf('Tc Dot size %.2f Gap Size %.2f',dotsize(ii),gapsize(jj)));
        colormap(jet);
        colorbar;
        
        figure;
        [DIFF, VL] = meshgrid(diff, vl);
        s = surf(DIFF, VL, tcOther(:,:)/sparams.ee*1E3);
        ylabel('V_L [V]');
        xlabel('V_L-V_R [V]');
        set(s,'edgecolor','none');
        view(2);
        xlim([min(diff),max(diff)]);
        title(sprintf('TcOther Dot size %.2f Gap Size %.2f',dotsize(ii),gapsize(jj)));
        colormap(jet);
        colorbar;
        drawnow;
        pp = pp + 1;
    end
end