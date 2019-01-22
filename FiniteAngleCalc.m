% Calculate mie scattering intensity 

% SEAFLOW scattering angles, -33.4:.1:-5.7, this NA = 0.55 using air for n1 [NA = n1 x sin(angle * pi/180)]
                           % -24.3:.1:-5.7, this NA = 0.55 using water for n1 
% [NA = n1 x sin(angle * pi/180)], so angle = asin(NA/n1) in radians, or angle = asin(NA/n1)* 180/pi in degree

%clear all;

tic	% start timer

%%%%  SEAFLOW PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%
m = 1.6;				% 1.017; 1.031; 1.045 refractive index ratio (air/water),1.6 beads, 1.3371 water
lam = 457e-9;			% SEAFLOW laser wavelength
k = 2*pi/lam;			% wavenumber
rmax=10;				% maximum radius in microns
r = exp(linspace(log(0.1),log(rmax),2000))*1e-6;	% particle radius
angle = [-33.4:.1:-5.7];   % SEAFLOW (NA = 0.55) scattering angles, 
x = k*r;				%
ang = angle*pi/180;		% put angle in radians for calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%  INFLUX PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%
% m = 1.6;				% 1.017; 1.031; 1.045 refractive index ratio (air/water),1.6 beads, 1.3371 water
% lam = 488e-9;			% SEAFLOW laser wavelength
% k = 2*pi/lam;			% wavenumber
% rmax=10;				% maximum radius in microns
% r = exp(linspace(log(0.1),log(rmax),2000))*1e-6;	% particle radius
% angle = [-18.30725:.1:-5.7];    % INFLUX  (NA = 0.42) scattering angles,                         
% x = k*r;				%
% ang = angle*pi/180;		% put angle in radians for calc
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% Ipara = parallel component of Intensity
% Iperp = perpendicular component of Intensity
[Ipara,Iperp] = Intensity(m, x, ang);	% call calculation subroutine

fprintf(1,'\nElapsed Time = %.0f seconds\n\n',toc)	% stop timer

% Average the angular intensities
avgIpara=sum(Ipara')/(length(ang));
avgIperp=sum(Iperp')/(length(ang));

% Plot results
figure;
semilogy(2*r*1e6,avgIpara,'r');grid on % multiply by 2 for diameter

%plot(angle,Iperp,'b',angle,Ipara,'r',...
%r*1e6,max(avgIpara)/rmax^2*(r*1e6).^2,'--k','LineWidth',2)
%legend('Perpendicular','Parallel','r^2 fit')
%plotaspect
ylabel('Relative Intensity');xlabel('diameter(µm)')
legend('Parallel','Perpendicular')
%title('cells n=1.063, 155-168')
%title('n=1.063, 6-33 degrees, NA = 0.55')



outPut=[2*r*1e6;avgIpara]; % output in diameter and micron
csvwrite('meidata-beadsb.csv',outPut);
% csvwrite('meidata-beadsINFLUX.csv',outPut);

