function [a,b,M]=Mmain(ftest4,hypotout)
    ev_mag=ftest4(:,4);
    mc=calc_McMaxCurvature(ev_mag);
    mag_sel = ev_mag(ev_mag >= mc+0.2);
    [b, bstd, a] =  calc_bmemag(mag_sel, 0.01);
    N=1;
    M=(log10(N)-a)/b;
end



