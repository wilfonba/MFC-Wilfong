#:def compute_capilary_tensor_x()

w1L = gL_x(j, k, l, 1)
w2L = gL_x(j, k, l, 2)
w3L = 0d0
if (p > 0) w3L = gL_x(j, k, l, 3)

w1R = gR_x(j + 1, k, l, 1)
w2R = gR_x(j + 1, k, l, 2)
w3R = 0d0
if (p > 0) w3R = gR_x(j + 1, k, l, 3)

normWL = gL_x(j, k, l, num_dims + 1)
normWR = gR_x(j + 1, k, l, num_dims + 1)

w1 = (w1L + w1R)/2d0
w2 = (w2L + w2R)/2d0
w3 = (w3L + w3R)/2d0
normW = (normWL + normWR)/2d0

if (normW > 1d-6) then
    Omega(1,1) = -sigma*(w2*w2 + w3*w3) / normW

    Omega(2,1) = sigma*w1*w2 / normW
    Omega(1,2) = Omega(2,1)

    Omega(2,2) = -sigma*(w1*w1 + w3*w3) / normW 
    
    if (p > 0) then

        Omega(3,1) = sigma*w1*w3 / normW
        Omega(1,3) = Omega(3,1)

        Omega(3,2) = sigma*w2*w2 / normW
        Omega(2,3) = Omega(3,2)

        Omega(3,3) = -sigma*(w1*w1 + w2*w2) / normW

    end if
else
    Omega(1,1) = 0d0
    Omega(1,2) = 0d0
    Omega(2,1) = 0d0
    Omega(2,2) = 0d0
    if (p > 0) then
        Omega(1,3) = 0d0
        Omega(3,1) = 0d0
        Omega(2,3) = 0d0
        Omega(3,2) = 0d0
        Omega(3,3) = 0d0
    end if
end if

#:enddef compute_capilary_tensor_x

#:def compute_capilary_tensor_y()

w1L = gL_y(k, j, l, 1)
w2L = gL_y(k, j, l, 2)
w3L = 0d0
if (p > 0) w3L = gL_y(k, j, l, 3)

w1R = gR_y(k + 1, j, l, 1)
w2R = gR_y(k + 1, j, l, 2)
w3R = 0d0
if (p > 0) w3R = gR_y(k + 1, j, l, 3)

normWL = gL_y(k, j, l, num_dims + 1)
normWR = gR_y(k + 1, j, l, num_dims + 1)

w1 = (w1L + w1R)/2d0
w2 = (w2L + w2R)/2d0
w3 = (w3L + w3R)/2d0
normW = (normWL + normWR)/2d0

if (normW > 1d-6) then
    Omega(1,1) = -sigma*(w2*w2 + w3*w3) / normW

    Omega(2,1) = sigma*w1*w2 / normW
    Omega(1,2) = Omega(2,1)

    Omega(2,2) = -sigma*(w1*w1 + w3*w3) / normW 
    
    if (p > 0) then

        Omega(3,1) = sigma*w1*w3 / normW
        Omega(1,3) = Omega(3,1)

        Omega(3,2) = sigma*w2*w2 / normW
        Omega(2,3) = Omega(3,2)

        Omega(3,3) = -sigma*(w1*w1 + w2*w2) / normW

    end if
else
    Omega(1,1) = 0d0
    Omega(1,2) = 0d0
    Omega(2,1) = 0d0
    Omega(2,2) = 0d0
    if (p > 0) then
        Omega(1,3) = 0d0
        Omega(3,1) = 0d0
        Omega(2,3) = 0d0
        Omega(3,2) = 0d0
        Omega(3,3) = 0d0
    end if
end if

#:enddef compute_capilary_tensor_y

#:def compute_capilary_tensor_z()

w1L = gL_z(l, j, k, 1)
w2L = gL_z(l, j, k, 2)
w3L = 0d0
if (p > 0) w3L = gL_z(l, j, k, 3)

w1R = gR_z(l + 1, j, k, 1)
w2R = gR_z(l + 1, j, k, 2)
w3R = 0d0
if (p > 0) w3R = gR_z(l + 1, j, k, 3)

normWL = gL_z(l, j, k, num_dims + 1)
normWR = gR_z(l + 1, j, k, num_dims + 1)

w1 = (w1L + w1R)/2d0
w2 = (w2L + w2R)/2d0
w3 = (w3L + w3R)/2d0
normW = (normWL + normWR)/2d0

if (normW > 1d-6) then
    Omega(1,1) = -sigma*(w2*w2 + w3*w3) / normW

    Omega(2,1) = sigma*w1*w2 / normW
    Omega(1,2) = Omega(2,1)

    Omega(2,2) = -sigma*(w1*w1 + w3*w3) / normW 
    
    if (p > 0) then

        Omega(3,1) = sigma*w1*w3 / normW
        Omega(1,3) = Omega(3,1)

        Omega(3,2) = sigma*w2*w2 / normW
        Omega(2,3) = Omega(3,2)

        Omega(3,3) = -sigma*(w1*w1 + w2*w2) / normW

    end if
else
    Omega(1,1) = 0d0
    Omega(1,2) = 0d0
    Omega(2,1) = 0d0
    Omega(2,2) = 0d0
    if (p > 0) then
        Omega(1,3) = 0d0
        Omega(3,1) = 0d0
        Omega(2,3) = 0d0
        Omega(3,2) = 0d0
        Omega(3,3) = 0d0
    end if
end if

#:enddef compute_capilary_tensor_z
