#:def compute_slope_limiter()

    if (muscl_lim == 1) then ! M
        phir = max(0d0, min(1d0, r))
    else if (muscl_lim == 2) then ! VanLeer
        if (r <= 0d0) phir = 0d0
        if (r >= 0d0) phir = min(2d0*r/(1 + r), 2d0/(1 + r))
    else if (muscl_lim == 3) then ! VanAlbada
        phir = max(0d0, (r + r**2d0)/(1d0 + r**2))
    else if (muscl_lim == 4) then ! Superbee
        phir = max(0d0, min(1d0, 2d0*r), min(2d0, r))
    end if

#:enddef

