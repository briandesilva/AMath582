%% Amath 582 Homework 3
% Brian de Silva
% 1242824

clear all; close all;

%% Task 1 - clean up the corrupted images of Derek

% % Global parameters (set specific filter parameters below)
generate_plots = true;                % Determines whether or not script makes plots
% filter_choice = 'threshold';           % Choice of filter: {'gaussian', 'shannon', 'threshold'}
% 
% 
% Acol = imread('derek1','jpg');      % Load color image
% Abw = imread('derek2','jpg');       % Load black & white image
% 
% % Convert uint8 arrays to double arrays and split the color image
% % into its RGB components
% Abw = double(Abw);
% AR = double(Acol(:,:,1));
% AG = double(Acol(:,:,2));
% AB = double(Acol(:,:,3));
% 
% % Take 2D FFT of each image/image component
% Abwt = fft2(Abw);
% ARt = fft2(AR);
% AGt = fft2(AG);
% ABt = fft2(AB);
% 
% 
% % Construct frequency vectors
% kxcenter = 181;                          % Specify center frequencies
% kycenter = 127;
% kx = 1:size(Abw,2);                     % Shifted frequencies
% ky = 1:size(Abw,1);
% [Kx, Ky] = meshgrid(kx,ky);
% 
% switch filter_choice
%     case 'gaussian'
%         % Contstruct Gaussian filter
%         xwidth = 0.0001;
%         ywidth = 0.0001;
%         filter = exp(-xwidth*(kxcenter - Kx).^2 - ywidth*(kycenter - Ky).^2);
%     case 'shannon'
%         % Construct Shannon (step function) filter
%         width = 50;
%         filter = zeros(size(Abwt));
%         filter((kycenter-width):(kycenter+width),(kxcenter-width):(kxcenter+width)) = ones(2*width+1,2*width+1);
%     case 'threshold'
%         % Filter out frequency components with too low magnitude
%         thresh = 0.6;
%         ARtf = ARt .* (log(abs(ARt)+1) > thresh * max(log(abs(ARt(:))+1)));
%         AGtf = AGt .* (log(abs(AGt)+1) > thresh * max(log(abs(AGt(:))+1)));
%         ABtf = ABt .* (log(abs(ABt)+1) > thresh * max(log(abs(ABt(:))+1)));
%         Abwtf = Abwt .* (log(abs(Abwt)+1) > thresh * max(log(abs(Abwt(:))+1)));
% end
% 
% % Apply filter (if not using thresholding)
% if exist('filter','var')
%     ARtf = ifftshift(filter.*fftshift(ARt));
%     AGtf = ifftshift(filter.*fftshift(AGt));
%     ABtf = ifftshift(filter.*fftshift(ABt));
%     Abwtf = ifftshift(filter.*fftshift(Abwt));
% end
% 
% % Plot Fourier components 
% if generate_plots
%     
%     % Color image
%     figure()
%     subplot(3,2,1)
%     pcolor(log(fftshift(abs(ARt)+1))), shading interp
%     title('log of magnitude of FFT for red component'), xlabel('Kx'), ylabel('Ky')
%     subplot(3,2,3)
%     pcolor(log(fftshift(abs(AGt))+1)), shading interp
%     title('log of magnitude of FFT for green component'), xlabel('Kx'), ylabel('Ky')
%     subplot(3,2,5)
%     pcolor(log(fftshift(abs(ABt))+1)), shading interp
%     title('log of magnitude of FFT for blue component'), xlabel('Kx'), ylabel('Ky')
%     subplot(3,2,2)
%     pcolor(log(fftshift(abs(ARtf))+1)), shading interp
%     title('log of magnitude of FFT for red component after filtering'), xlabel('Kx'), ylabel('Ky')
%     subplot(3,2,4)
%     pcolor(log(fftshift(abs(AGtf))+1)), shading interp
%     title('log of magnitude of FFT for green component after filtering'), xlabel('Kx'), ylabel('Ky')
%     subplot(3,2,6)
%     pcolor(log(fftshift(abs(ABtf))+1)), shading interp
%     title('log of magnitude of FFT for blue component after filtering'), xlabel('Kx'), ylabel('Ky')
% 
%     % Black and white image
%     figure()
%     subplot(2,1,1)
%     pcolor(log(fftshift(abs(Abwt))+1)), shading interp
%     title('log of magnitude of FFT for bw image'), xlabel('Kx'), ylabel('Ky')
%     subplot(2,1,2)
%     pcolor(log(fftshift(abs(Abwtf))+1)), shading interp
%     title('log of magnitude of FFT for bw image after filtering'), xlabel('Kx'), ylabel('Ky')
% end
% 
% 
% % Invert transform
% 
% % Overwrite original image to save space
% % Acol(:,:,1) = uint8(ifft2(ARtf));
% % Acol(:,:,2) = uint8(ifft2(AGtf));
% % Acol(:,:,3) = uint8(ifft2(ABtf));
% % Abwf = uint8(ifft2(Abwtf));
% 
% Acolf = uint8(zeros(size(Acol)));
% Acolf(:,:,1) = uint8(ifft2(ARtf));
% Acolf(:,:,2) = uint8(ifft2(AGtf));
% Acolf(:,:,3) = uint8(ifft2(ABtf));
% Abwf = uint8(ifft2(Abwtf));
% 
% % Plot original and filtered images
% if generate_plots
%     % Color images
%     figure()
%     subplot(2,1,1)
%     image(Acol); set(gca,'Xtick',[],'Ytick',[])
%     title('Sabotaged color photo of Derek')
%     
%     subplot(2,1,2)
%     image(Acolf); set(gca,'Xtick',[],'Ytick',[])
%     title('Filtered color photo of Derek')
% 
%     % Black and white images
%     figure()
%     subplot(2,1,1);
%     imshow(uint8(Abw))
%     title('Sabotaged black and white photo of Derek')
% %     title(findobj('Parent',hFig,'Type','axes'),'Sabotaged black and white photo of Derek')
%     subplot(2,1,2);
%     imshow(Abwf)
%     title('Filtered black and white photo of Derek')
% %     title(findobj('Parent',hFig,'Type','axes'),'Filtered black and white photo of Derek');
% end



%% Task 2: Remove Derek's rash

Ac = imread('derek3','jpg');          % Load color image
AR = Ac(:,:,1);                              % Get red component of color image
AG = Ac(:,:,2);                              % Get green component of color image
AB = Ac(:,:,3);                              % Get blue component of color image
Abw = imread('derek4','jpg');       % Load black and white image
[nx, ny] = size(Abw);



% Indices of rash: (140:160, 160:180)

% Construct discretization of the Laplacian
x = linspace(0,1,nx);
y = linspace(0,1,ny);
dx = x(2) - x(1);
dy = y(2) - y(1);
onex = ones(nx,1);
oney = ones(ny,1);
Dx = spdiags([onex, -2*onex, onex],[-1 0 1],nx,nx) / (dx^2);
Ix = speye(nx);
Dy = spdiags([oney, -2*oney, oney],[-1 0 1],ny,ny) / (dy^2);
Iy = speye(ny);
L = kron(Iy,Dx) + kron(Dy, Ix);        % Discretization of the Laplacian
D = 0.0005;                                     % Diffusion coefficient inside rash domain

% Create matrix to enforce that no diffusion happens outside of the rash
B = zeros(size(Abw));
B(140:160,160:180) = 1;

% Define function encapsulating rhs of the ode to be solved
odefun = @(t,u) (D*L*u) .* B(:);

% Solve the ODE at different times
tpts = [0, 0.015, 0.03, 0.045];                                          % Points in time when we want to view solution
[~, usolnR] = ode45(odefun, tpts, double(AR(:)));               %  Solve the ODE for the red component
[~, usolnG] = ode45(odefun, tpts, double(AG(:)));              %  Solve the ODE for the green component
[~, usolnB] = ode45(odefun, tpts, double(AB(:)));              %  Solve the ODE for the blue component
[~, usolnbw] = ode45(odefun, tpts, double(Abw(:)));         % Solve the ODE for the black and white image

% Plot results
if generate_plots
    % Color
%     figure()
    for i=1:4
%         subplot(2,2,i)
        figure()
        image(uint8(reshape(cat(3,usolnR(i,:),usolnG(i,:),usolnB(i,:)),nx,ny,3)));
        set(gca,'Xtick',[],'Ytick',[]);
        title(['t=' num2str(tpts(i))]);
    end
    
    % Black and white
%     figure()
    for i=1:4
%         subplot(2,2,i)
        figure()
        imshow(uint8(reshape(usolnbw(i,:),nx,ny)));
        title(['t=' num2str(tpts(i))]);
    end

end

















