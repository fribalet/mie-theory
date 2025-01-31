% Calculate mie scattering intensity from bubble in water
%
clear all;
tic	% start timer
%%%%  PARAMETERS  %%%%%%%%%%%%%%%%%%%%%%
m = 1.035/1.341;				% refractive index ratio (air/water), 1.6033 beads 1.035 cells, 1.33 water
lam = 457e-9;			% laser wavelength
k = 2*pi/lam;			% wavenumber
%rmax=100;				% maximum radius in microns
r = 6*0.25e-6;			% radius
angle = [0:.01:6];	% scattering angles
x = k*r;				% 
ang = angle.*pi/180;		% put angle in radians for calc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ipara = parallel component of Intensity
% Iperp = perpendicular component of Intensity
[Ipara,Iperp] = Intensity(m, x, ang);	% call calculation subroutine

% calculate differentiation of intensity
Num_point=size(angle,2);
Ipara_d=zeros(1,Num_point-1);
Iperp_d=zeros(1,Num_point-1);
for nn=2:Num_point
   Ipara_d(nn-1)=Ipara(nn)-Ipara(nn-1);
   Iperp_d(nn-1)=Iperp(nn)-Iperp(nn-1);
end

fprintf(1,'\nElapsed Time = %.0f seconds\n\n',toc)	% stop timer

% Average the angular intensities
%avgIpara=sum(Ipara');
%avgIperp=sum(Iperp');

% Plot results
% figure(1);
% plot(angle,Iperp,'b',angle,Ipara,'r')
% legend('Perpendicular','Parallel')
% %plotaspect
% ylabel('Relative Intensity');xlabel('angle (degree)')

figure(2);
plot(angle,log10(Iperp),'b',angle,log10(Ipara),'r')
legend('Perpendicular','Parallel')
%plotaspect
ylabel('Relative Intensity');xlabel('angle (degree)')

% figure(3);
% angle2=angle;
% angle2(Num_point)=[];
% plot(angle2,Iperp_d,'b',angle2,Ipara_d,'r')
% legend('Perpendicular','Parallel')
% %plotaspect
% ylabel('Relative differentiated Intensity');xlabel('angle (degree)')


figure
polar(ang,Ipara)
title('6 um bead scatter, 60-180 degrees')
