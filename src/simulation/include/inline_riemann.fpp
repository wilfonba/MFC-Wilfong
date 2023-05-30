#:def arithmetic_avg()
    rho_avg = 5d-1*(rho_L + rho_R)
    vel_avg_rms = 0d0
    !$acc loop seq
    do i = 1, num_dims
        vel_avg_rms = vel_avg_rms + (5d-1*(vel_L(i) + vel_R(i)))**2d0
    end do

    H_avg = 5d-1*(H_L + H_R)
    gamma_avg = 5d-1*(gamma_L + gamma_R)

#:enddef arithmetic_avg


#:def roe_avg()
    rho_avg = sqrt(rho_L*rho_R)
    vel_avg_rms = 0d0
    !$acc loop seq
    do i = 1, num_dims
        vel_avg_rms = vel_avg_rms + (sqrt(rho_L)*vel_L(i) + sqrt(rho_R)*vel_R(i))**2d0/ &
                      (sqrt(rho_L) + sqrt(rho_R))**2d0
    end do

    H_avg = (sqrt(rho_L)*H_L + sqrt(rho_R)*H_R)/ &
            (sqrt(rho_L) + sqrt(rho_R))

    gamma_avg = (sqrt(rho_L)*gamma_L + sqrt(rho_R)*gamma_R)/ &
                (sqrt(rho_L) + sqrt(rho_R))

    rho_avg = sqrt(rho_L*rho_R)
    vel_avg_rms = (sqrt(rho_L)*vel_L(1) + sqrt(rho_R)*vel_R(1))**2d0/ &
                  (sqrt(rho_L) + sqrt(rho_R))**2d0

#:enddef roe_avg

#:def lo_arithmetic_avg()

    lo_rho_avg = 5d-1*(lo_rho_L + lo_rho_R)
    lo_vel_avg_rms = 0d0
    !$acc loop seq
    do i = 1, num_dims
        lo_vel_avg_rms = lo_vel_avg_rms + (5d-1*(lo_vel_L(i) + lo_vel_R(i)))**2d0
    end do

    lo_H_avg = 5d-1*(lo_H_L + lo_H_R)
    lo_gamma_avg = 5d-1*(lo_gamma_L + lo_gamma_R)

#:enddef lo_arithmetic_avg


#:def lo_roe_avg()
    lo_rho_avg = sqrt(lo_rho_L*lo_rho_R)
    lo_vel_avg_rms = 0d0
    !$acc loop seq
    do i = 1, num_dims
        lo_vel_avg_rms = lo_vel_avg_rms + (sqrt(lo_rho_L)*lo_vel_L(i) + sqrt(lo_rho_R)*lo_vel_R(i))**2d0/ &
                      (sqrt(lo_rho_L) + sqrt(lo_rho_R))**2d0
    end do

    lo_H_avg = (sqrt(lo_rho_L)*lo_H_L + sqrt(lo_rho_R)*lo_H_R)/ &
            (sqrt(lo_rho_L) + sqrt(lo_rho_R))

    lo_gamma_avg = (sqrt(lo_rho_L)*lo_gamma_L + sqrt(lo_rho_R)*lo_gamma_R)/ &
                (sqrt(lo_rho_L) + sqrt(lo_rho_R))

    lo_rho_avg = sqrt(lo_rho_L*lo_rho_R)
    lo_vel_avg_rms = (sqrt(lo_rho_L)*lo_vel_L(1) + sqrt(lo_rho_R)*lo_vel_R(1))**2d0/ &
                  (sqrt(lo_rho_L) + sqrt(lo_rho_R))**2d0

#:enddef lo_roe_avg


#:def compute_average_state()

if (avg_state == 1) then
    @:roe_avg()
end if

if (avg_state == 2) then
    @:arithmetic_avg()
end if

#:enddef compute_average_state


#:def lo_compute_average_state()

if (avg_state == 1) then
    @:lo_roe_avg()
end if

if (avg_state == 2) then
    @:lo_arithmetic_avg()
end if

#:enddef lo_compute_average_state

