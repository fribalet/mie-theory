% CROSSSECTION	Scattering & Extinction Cross Sections%		[Cext, Cscat] = crosssection(m,a,wv) returns the extinction%		and scattering cross sections for a sphere, radius a,%		refractive index relative to medium m at wavelength wv%%		If any of the parameters are vectors, return will a vector.% Written by and copyright%		Dave Barnett%		Optical Engineering Group%		Dept. of Electrical and Electronic Engineering%		Loughborough University%		4th May 1999function [Cext, Cscat] = crosssection(m,a,wv)wv = wv(:).';a = a(:).';m = m(:).';x = 2*pi*a./wv;if length(x)==1   x = x*ones(size(m));endif length(m)==1   m = m*ones(size(x));endnc = ceil(max(x)+4.05*(max(x)^(1/3))+2);n=(1:nc).';W = warning;warning off[a,b] = ScatCoef(m,x,nc);% Check for invalid (NaN) results due to too many terms in% relatively small particles.invalid = find(any(isnan([a;b])));while ~isempty(invalid)   a(:,invalid) = 0;   b(:,invalid) = 0;   nc2 = ceil(max(x(invalid))+4.05*(max(x(invalid))^(1/3))+2);   [A,B] = ScatCoef(m(invalid),x(invalid),nc2);   a(1:nc2,invalid) = A;   b(1:nc2,invalid) = B;   invalid = find(any(isnan([a;b])));   % remove invalidity of zero m or x   % these _should_ return NaN!   if length(x)>=max(invalid)      invalid = invalid(x(invalid)~=0);   else      if x==0         invalid = [];      end   end   if length(m)>=max(invalid)      invalid = invalid(m(invalid)~=0);   else      if m==0         invalid = [];      end   endendwarning(W);Cext = ((wv.^2)/(2*pi)).*((2*n+1)'*real(a+b));Cscat = ((wv.^2)/(2*pi)).*((2*n+1)'*(abs(a).^2+abs(b).^2));