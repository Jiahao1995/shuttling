% Load all the parameters for the simulation
clear sparams xx vv;
shuttleParameterFile;

fprintf(1,'Loading potentials...\n');
[xx,zz,vv] = loadPotentials(sparams);

sparams.nxGrid = length(xx);
sparams.dx = xx(2) - xx(1);
sparams.dp = 2*pi*sparams.hbar/(sparams.dx*sparams.nxGrid);
pp = ((-sparams.nxGrid/2):1:(sparams.nxGrid/2 - 1))*sparams.dp;

% Find which index corresponds to where the 2DEG should be
[~,index] = min(abs(zz - (-0.5*1E-9)));
sparams.twoDEGindZ = index;
sparams.twoDEGPots = squeeze(vv(:,sparams.twoDEGindZ,:));

% Check that the potentials were loaded correctly
[nPots,~,~] = size(vv);
debugHere = 1;
if debugHere
    fig = figure;
    
    yminlim = min(min(sparams.twoDEGPots));
    ymaxlim = max(max(sparams.twoDEGPots));
    
    for ii = 1:nPots
        clf;
        oneDPotSlice = sparams.twoDEGPots(ii,:);
        plot(xx,oneDPotSlice);
        ylim([yminlim, ymaxlim]);
        pause(3/100);
    end
    delete(fig);
    
    % Plot 2D version of 1st potential to make sure it was imported
    % correctly
    [XX,ZZ] = meshgrid(xx,zz);
    figure;
    s = surf(XX,ZZ,squeeze(vv(1,:,:)));
    set(s,'edgecolor','none');
    view(2);
end

fprintf(1,'Getting initial wavefunction...\n');
% Solve the 1D SE for the initial potential well to get what our ground
% state should look like
[sparams.rho0, ~] = solve1DSingleElectronSE(sparams,1,xx,sparams.twoDEGPots(1,:)); 

% Check that the intial state makes sense
debugHere = 1;
if debugHere
    figure;
    hold on;
    plot(xx,sparams.twoDEGPots(1,:)/sparams.ee);
    plot(xx,abs(sparams.rho0).^2 + min(sparams.twoDEGPots(1,:)/sparams.ee));
    title('Initial conditions');
end
%%
% Using ref https://arxiv.org/pdf/1306.3247.pdf we now find the time
% evolution operator U(t + dT,t) to evolve our initial wavefunction to the
% next wavefunction.  That newly found wavefunction will act as our initial
% state for the next time frame.  This method uses the split ooperator
% approach
    
% Make the KE operator since it is the same every time it is applied
K = exp(-1i*sparams.dt/2*(pp.^2)/(2*sparams.me*sparams.hbar));

% Make the fidelity array
maxTime = max(sparams.totalTime);
maxLength = length(0:sparams.dt:maxTime);
sparams.fidelity = zeros(length(sparams.totalTime),floor(maxLength/sparams.updateFidelity));

for jj = 1:length(sparams.totalTime)
    tic;
    
    % Now, we want to associate each potential simulation we have with a time
    % value (i.e. when in the simulation should that potential appear)
    [nPots,~] = size(vv);
%     tPots = linspace(0,sparams.totalTime(jj),nPots);
    NP = 16;
    tPots = linspace(0,sparams.totalTime(jj)*1/16,nPots/4 - NP);
    tPotsTemp = linspace(sparams.totalTime(jj)*1/16,sparams.totalTime(jj)*7/16,2*NP + 1);
    tPots = [tPots, tPotsTemp(2:end)];
    tPotsTemp = linspace(sparams.totalTime(jj)*7/16,sparams.totalTime(jj)*9/16,nPots/2 - 2*NP + 1);
    tPots = [tPots, tPotsTemp(2:end)];
    tPotsTemp = linspace(sparams.totalTime(jj)*9/16,sparams.totalTime(jj)*15/16,2*NP + 1);
    tPots = [tPots, tPotsTemp(2:end)];
    tPotsTemp = linspace(sparams.totalTime(jj)*15/16,sparams.totalTime(jj)*16/16,nPots/4 - NP + 1);
    tPots = [tPots, tPotsTemp(2:end)];
    
    % Get number of time steps
    tTime = 0:sparams.dt:sparams.totalTime(jj);
    nTime = length(tTime);

    saveFigureIndices = 1:round(nTime/sparams.nFigureFrames):nTime;
    if saveFigureIndices(end) ~= nTime
        saveFigureIndices = [saveFigureIndices nTime];
    end
    
    fprintf(1,'Running shuttling simulation for %E (%d/%d)...\n',sparams.totalTime(jj),jj,length(sparams.totalTime));

    % Make the waitbar to show run time
    str1 = sprintf('Current Time Index: %d/%d',0,nTime);
    str2 = sprintf('Performing shuttling simulation for %E...',sparams.totalTime(jj));
    h = waitbar(0,str1,'Name',str2,'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    movegui(h,'northwest');

    currPsi = sparams.rho0';
    currFig = figure;
    nn = 1; % Used to index fidelity array
    mm = 1; % Used to index gif creation
    ll = 0; % Used to know where in time domain to interpolate our potentials
    kk = 0; % Used to index which interpolated potential we are on

    for ii = 1:nTime
        kk = kk + 1;
        
        %Check for cancel button click
        if getappdata(h,'canceling')
            flag = 1;
            break;
        end

        % Update waitbar every N frames
        if mod(ii,sparams.updateWaitbar) == 0
            waitbar(ii/nTime, h, sprintf('Current Time Index: %d/%d',ii,nTime));
        end

        % Get updated set of interpolated potentials if needed
        if strcmp(sparams.interpType,'linear')
            if mod(ii,sparams.updateInterpPot) == 0 || ii == 1
                kk = 1; % Reset counter
                startInterpInd = ll*sparams.updateInterpPot;
                if ii == 1
                    startInterpInd = 1;
                end
                ll = ll + 1;
                endInterpInd = ll*sparams.updateInterpPot - 1;
                if endInterpInd > nTime
                    endInterpInd = nTime;
                end
                vvInterp = interp1(tPots,sparams.twoDEGPots,tTime(startInterpInd:endInterpInd));
            end
        end
        V = exp(-1i*sparams.dt*vvInterp(kk,:)/sparams.hbar);

        % Every frame, we want to interpolate what the potential should look
        % like at that particular time
    %     if strcmp(sparams.interpType,'linear')
    %         [time,ind] = min(abs(sparams.tPots - sparams.tTime(ii)));
    %         % See if our current time is less than or greater than or equal to
    %         % the closest index
    %         if sparams.tTime(ii) == sparams.tPots(ind);
    %             currPotential = vv(ind,:);
    %         elseif sparams.tTime(ii) < sparams.tPots(ind);
    %             currPotential = vv(ind-1,:) + (sparams.tTime(ii) - sparams.tPots(ind-1))*...
    %                 ((vv(ind,:) - vv(ind-1,:))/(sparams.tPots(ind) - sparams.tPots(ind-1)));
    %         elseif sparams.tTime(ii) > sparams.tPots(ind);
    %             currPotential = vv(ind,:) + (sparams.tTime(ii) - sparams.tPots(ind))*...
    %                 ((vv(ind+1,:) - vv(ind,:))/(sparams.tPots(ind+1) - sparams.tPots(ind)));
    %         end
    %     end
    %     V = exp(-1i*sparams.dt*currPotential/sparams.hbar);

        % Convert from position to momentum space
        currPsip = fftshift(fft(fftshift(currPsi)));
        % Apply the KE operator for dt/2 
        currPsip = K.*currPsip;
        % Convert from momentum to position space
        currPsix = fftshift(ifft(fftshift(currPsip)));
        % Apply the PE operator for dt
        currPsix = V.*currPsix;
        % Convert from position to momentum space
        currPsip = fftshift(fft(fftshift(currPsix)));
        % Apply the KE operator for dt/2
        currPsip = K.*currPsip;
        % Convert from momentum to position space
        currPsi = fftshift(ifft(fftshift(currPsip)));

        % Update figure every N frames and save to gif
        if mod(ii,sparams.updateFigure) == 0
            [currRho0, ~] = solve1DSingleElectronSE(sparams,1,xx,vvInterp(kk,:));
            clf;
            hold on;
            plot(xx,vvInterp(kk,:)/sparams.ee);
            plot(xx,2*abs(currPsi).^2 + min(sparams.twoDEGPots(1,:)/sparams.ee));
            plot(xx,2*abs(currRho0).^2 + min(sparams.twoDEGPots(1,:)/sparams.ee));
            title(['Shuttling Simulation ' num2str(sparams.totalTime(jj)) '[s]'],'interpreter','latex','fontsize',12);
            xlabel('Position [m]','interpreter','latex','fontsize',12);
            ylabel('Energy [eV]','interpreter','latex','fontsize',12);
%             ylim([min(min(vv)),max(max(vv))]);
        end

        % Update figure and save to gif according to figure frames
        % parameter
        if any(saveFigureIndices == ii)
            [currRho0, ~] = solve1DSingleElectronSE(sparams,1,xx,vvInterp(kk,:));
            clf;
            hold on;
            plot(xx,vvInterp(kk,:)/sparams.ee);
            plot(xx,2*abs(currPsi).^2 + min(sparams.twoDEGPots(1,:)/sparams.ee));
            plot(xx,2*abs(currRho0).^2 + min(sparams.twoDEGPots(1,:)/sparams.ee));
            title(['Shuttling Simulation ' num2str(sparams.totalTime(jj)) '[s]'],'interpreter','latex','fontsize',12);
            xlabel('Position [m]','interpreter','latex','fontsize',12);
            ylabel('Energy [eV]','interpreter','latex','fontsize',12);
%             ylim([min(min(vv)),max(max(vv))]);;
            currFrame = getframe(currFig);
            currIm{mm} = frame2im(currFrame);
            [A,map] = rgb2ind(currIm{mm},256);
            if mm == 1
                imwrite(A,map,['shuttle' num2str(sparams.totalTime(jj)) '.gif'],'gif','LoopCount',Inf,'DelayTime',0);
            else
                imwrite(A,map,['shuttle' num2str(sparams.totalTime(jj)) '.gif'],'gif','WriteMode','append','DelayTime',0);
            end
            mm = mm + 1;
        end
        
        % Calculate fidelity with current ground state every N frames
        if mod(ii,sparams.updateFidelity) == 0
            % Need to get the ground state of the current potential
    %         [currRho0, ~] = solve1DSingleElectronSE(sparams,1,xx,currPotential);
            [currRho0, ~] = solve1DSingleElectronSE(sparams,1,xx,vvInterp(kk,:));
            sparams.fidelity(jj,nn) = abs(sum(currRho0'.*currPsi))^2;
            nn = nn + 1;
        end
    end
    % Close simullation figure
    close(currFig);
    % Close waitbar
    delete(h);
    % Delete currIm
    clearvars currIm
    toc;
end


%%
% fids = sparams.fidelity;
fids = sparams.fidelity;
fids(find(fids==0)) = NaN;
[rows,cols] = size(fids);
highTime = max(sparams.totalTime);
fidTimeIndices = sparams.updateFidelity*sparams.dt:sparams.updateFidelity*sparams.dt:highTime;
[TIndex,TTime] = meshgrid(fidTimeIndices,[0,sparams.totalTime,2*max(sparams.totalTime)]);
fidelTemp = zeros(rows+2,cols);
fidelTemp(2:(rows+1),:) = fids;


figure;
% fidelTemp(fidelTemp == 0) = NaN;
s = surf(TIndex,TTime,fidelTemp);
set(s,'edgecolor','none');
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('Time step $t_j$ [s]','interpreter','latex','fontsize',15);
ylabel('Total Shuttling Simulated Time [s]','interpreter','latex','fontsize',15);
xlim([min(min(TIndex)),max(max(TIndex))]);
ylim([0,2*max(sparams.totalTime)]);
title('Fidelity: $$|\langle\Psi_0(t_j)|\Psi_{\rm sim}(t_j)\rangle|^2$$','interpreter','latex','fontsize',15);
view(2);
colormap(jet);
colorbar;


