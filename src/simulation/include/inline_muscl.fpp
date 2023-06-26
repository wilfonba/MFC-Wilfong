#:def compute_slope_limiter()

    if (muscl_lim == 1) then ! Van Albada
        phir = (r + r ** 2)/(1 + r ** 2)
    else if (muscl_lim == 2) then ! Van Leer
        phir = max(0d0, min((3d0/2d0)*r, (1d0 + r)/2d0, 3d0/2d0))
    else if (muscl_lim == 3) then ! Superbee
        phir = max(0d0, min(min(2*r,1d0), min(r,2d0)))
    end if

#:enddef

