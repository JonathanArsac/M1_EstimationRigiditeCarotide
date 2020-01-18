function [f,pixel_size_ax,pixel_size_lat]=compute_parameters_for_sta(x,z,lambdax,decim_factor)
pixel_size_ax=decim_factor*(max(z)-min(z))/length(z)*1e3; %dimension in mm
pixel_size_lat=decim_factor*(max(x)-min(x))/length(x)*1e3;

fc=5e6;
c=1540;
lambdaz=c/(2*fc);
f_ax=1e-3*pixel_size_ax/lambdaz;

f_lat=2e-3*pixel_size_lat/lambdax; %sta specific

f=[f_ax,f_lat];