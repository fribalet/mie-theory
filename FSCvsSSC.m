% Calculate mie scattering intensity from bubble in water
%
clear all;
tic	% start timer
%%%%  PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%
m = 1.035/1.3371;				% refractive index ratio (air/water)1.6033 for beads, 1.035 for plankton
lam = 488e-9;			% laser wavelength
k = 2*pi/lam;			% wavenumber
rmax=10;				% maximum radius in microns
r = [1:.25:rmax]'*1e-6;	% bubble radius
angle = [-8:.1:8]; %angle = 80;	% scattering angles
x = k*r;				% 
ang = angle*pi/180;		%put angle in radians for calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ipara = parallel component of Intensity
% Iperp = perpendicular component of Intensity
[Ipara,Iperp] = Intensity(m, x, ang);	% call calculation subroutine

fprintf(1,'\nElapsed Time = %.0f seconds\n\n',toc)	% stop timer

% Average the angular intensities
avgIpara=sum(Ipara');
avgIperp=sum(Iperp');

% Plot results
figure;
loglog(r*1e6,avgIperp,'b',r*1e6,avgIpara,'r');grid on
%plot(angle,Iperp,'b',angle,Ipara,'r',...
%r*1e6,max(avgIpara)/rmax^2*(r*1e6).^2,'--k','LineWidth',2)
%legend('Perpendicular','Parallel','r^2 fit')
%plotaspect
ylabel('Relative Intensity');xlabel('radius (µm)')
legend('Perpendicular','Parallel')
title('FSC')
outPut=[r'*1e6;avgIperp];

%fid=fopen('C:\users\Jing\Matlab program\ContrastAgent\MieIR1.txt', 'w');
%fprintf(fid,'%12.8f     %12.8f\n',outPut);
%fclose(fid);




%Second optical path


angleSSC = [82:0.1:98]; %angle = 80;	% scattering angles
xSSC = k*r;				 
angSSC = angleSSC*pi/180;		

[IparaSSC,IperpSSC] = Intensity(m, xSSC, angSSC);	% call calculation subroutine

fprintf(1,'\nElapsed Time = %.0f seconds\n\n',toc)	% stop timer

avgIparaSSC=sum(IparaSSC');
avgIperpSSC=sum(IperpSSC');

figure;
loglog(r*1e6,avgIperpSSC,'b',r*1e6,avgIparaSSC,'r');grid on
ylabel('Relative Intensity');xlabel('radius (µm)')
legend('Perpendicular','Parallel')
title('SSC')

outPut=[r'*1e6;avgIperp];

figure;
loglog(avgIperp,avgIperpSSC,'r.')
ylabel('Perp SSC')
xlabel('Perp FSC')
title('Perpindicular FSC vs. SSC')

figure
loglog(avgIpara,avgIparaSSC,'r.')
ylabel('Para SSC')
xlabel('Para FSC')
title('Parallel FSC vs. SSC')

