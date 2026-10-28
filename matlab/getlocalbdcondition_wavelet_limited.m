function localbd = getlocalbdcondition_wavelet_limited(nx, ny, type, ncols)
%GETLOCALBDCONDITION_WAVELET_LIMITED Boundary wavelets needed by this study.
%   getlocalbdcondition_select builds 64 columns.  For the Darcy wavelet
%   sweep here we only need columns 1:4, 1:8, or 1:16.  Avoid constructing
%   the higher columns because small oversampling patches need not be
%   divisible by 8 or 16.

localbd = zeros(2 * nx + 2 * ny, ncols);
if ncols <= 0
    return;
end
if ncols > 16
    error('This limited helper supports at most 16 boundary columns.');
end

if type == 1
    if ncols >= 1
        localbd(1:nx + 1, 1) = 1:-1 / nx:0;
        localbd(2 * nx + ny + 1:end, 1) = 0:1 / ny:1 - 1 / ny;
    end
    if ncols >= 2
        localbd(1:nx + 1, 2) = 0:1 / nx:1;
        localbd(nx + 1:nx + ny + 1, 2) = 1:-1 / ny:0;
    end
    if ncols >= 3
        localbd(nx + 1:nx + ny + 1, 3) = 0:1 / ny:1;
        localbd(nx + ny + 1:2 * nx + ny + 1, 3) = 1:-1 / nx:0;
    end
    if ncols >= 4
        localbd(nx + ny + 1:2 * nx + ny + 1, 4) = 0:1 / nx:1;
        localbd(2 * nx + ny + 1:end, 4) = 1:-1 / ny:1 / ny;
    end
    if ncols >= 5
        a = 0;
        n1 = floor(nx / 2);
        n2 = nx - n1;
        localbd(a + 1:a + n1 + 1, 5) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + nx + 1, 5) = 1:-1 / n2:0;
    end
    if ncols >= 6
        a = nx;
        n1 = floor(ny / 2);
        n2 = ny - n1;
        localbd(a + 1:a + n1 + 1, 6) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + ny + 1, 6) = 1:-1 / n2:0;
    end
    if ncols >= 7
        a = nx + ny;
        n1 = floor(nx / 2);
        n2 = nx - n1;
        localbd(a + 1:a + n1 + 1, 7) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + nx + 1, 7) = 1:-1 / n2:0;
    end
    if ncols >= 8
        a = 2 * nx + ny;
        n1 = floor(ny / 2);
        n2 = ny - n1;
        localbd(a + 1:a + n1 + 1, 8) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + ny, 8) = 1:-1 / n2:1 / n2;
    end
    if ncols >= 9
        a = 0;
        nt = floor(nx / 2);
        n1 = floor(nx / 4);
        n2 = nt - n1;
        localbd(a + 1:a + n1 + 1, 9) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + nt + 1, 9) = 1:-1 / n2:0;
    end
    if ncols >= 10
        a = floor(nx / 2);
        nt = nx - floor(nx / 2);
        n1 = floor(nx / 4);
        n2 = nt - n1;
        localbd(a + 1:a + n1 + 1, 10) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + nt + 1, 10) = 1:-1 / n2:0;
    end
    if ncols >= 11
        a = nx;
        nt = floor(ny / 2);
        n1 = floor(ny / 4);
        n2 = nt - n1;
        localbd(a + 1:a + n1 + 1, 11) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + nt + 1, 11) = 1:-1 / n2:0;
    end
    if ncols >= 12
        a = nx + floor(ny / 2);
        nt = nx - floor(nx / 2);
        n1 = floor(nx / 4);
        n2 = nt - n1;
        localbd(a + 1:a + n1 + 1, 12) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + nt + 1, 12) = 1:-1 / n2:0;
    end
    if ncols >= 13
        a = nx + ny;
        nt = floor(nx / 2);
        n1 = floor(nx / 4);
        n2 = nt - n1;
        localbd(a + 1:a + n1 + 1, 13) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + nt + 1, 13) = 1:-1 / n2:0;
    end
    if ncols >= 14
        a = nx + ny + floor(nx / 2);
        nt = floor(nx / 2);
        n1 = floor(nx / 4);
        n2 = nt - n1;
        localbd(a + 1:a + n1 + 1, 14) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + nt + 1, 14) = 1:-1 / n2:0;
    end
    if ncols >= 15
        a = nx + ny + nx;
        nt = floor(ny / 2);
        n1 = floor(ny / 4);
        n2 = nt - n1;
        localbd(a + 1:a + n1 + 1, 15) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + nt + 1, 15) = 1:-1 / n2:0;
    end
    if ncols >= 16
        a = nx + ny + nx + floor(ny / 2);
        nt = floor(ny / 2);
        n1 = floor(ny / 4);
        n2 = nt - n1;
        localbd(a + 1:a + n1 + 1, 16) = 0:1 / n1:1;
        localbd(a + n1 + 1:a + nt, 16) = 1:-1 / n2:1 / n2;
    end
else
    if ncols >= 1
        localbd(1:nx, 1) = 1;
    end
    if ncols >= 2
        localbd(nx + 1:nx + ny, 2) = 1;
    end
    if ncols >= 3
        localbd(nx + ny + 1:2 * nx + ny, 3) = 1;
    end
    if ncols >= 4
        localbd(2 * nx + ny + 1:end, 4) = 1;
    end
    if ncols >= 5
        localbd(1:nx / 2, 5) = 1;
        localbd(nx / 2 + 1:nx, 5) = -1;
    end
    if ncols >= 6
        localbd(1 + nx:ny / 2 + nx, 6) = 1;
        localbd(ny / 2 + 1 + nx:ny + nx, 6) = -1;
    end
    if ncols >= 7
        localbd(1 + nx + ny:nx / 2 + nx + ny, 7) = 1;
        localbd(nx / 2 + 1 + nx + ny:nx + nx + ny, 7) = -1;
    end
    if ncols >= 8
        localbd(1 + 2 * nx + ny:ny / 2 + 2 * nx + ny, 8) = 1;
        localbd(ny / 2 + 1 + 2 * nx + ny:ny + 2 * nx + ny, 8) = -1;
    end
    if ncols >= 9
        localbd(1:nx / 4, 9) = 1;
        localbd(nx / 4 + 1:nx / 2, 9) = -1;
    end
    if ncols >= 10
        localbd(1 + nx:ny / 4 + nx, 10) = 1;
        localbd(ny / 4 + 1 + nx:ny / 2 + nx, 10) = -1;
    end
    if ncols >= 11
        localbd(1 + nx + ny:nx / 4 + nx + ny, 11) = 1;
        localbd(nx / 4 + 1 + nx + ny:nx / 2 + nx + ny, 11) = -1;
    end
    if ncols >= 12
        localbd(1 + nx * 2 + ny:ny / 4 + nx * 2 + ny, 12) = 1;
        localbd(nx / 4 + 1 + nx * 2 + ny:ny / 2 + nx * 2 + ny, 12) = -1;
    end
    if ncols >= 13
        localbd(1 + nx / 2:3 / 4 * nx, 13) = 1;
        localbd(3 * nx / 4 + 1:nx, 13) = -1;
    end
    if ncols >= 14
        localbd(1 + ny / 2 + nx:3 * ny / 4 + nx, 14) = 1;
        localbd(3 * ny / 4 + 1 + nx:ny + nx, 14) = -1;
    end
    if ncols >= 15
        localbd(1 + nx / 2 + nx + ny:3 * nx / 4 + nx + ny, 15) = 1;
        localbd(3 * nx / 4 + 1 + nx + ny:nx + nx + ny, 15) = -1;
    end
    if ncols >= 16
        localbd(1 + nx / 2 + nx * 2 + ny:3 * ny / 4 + nx * 2 + ny, 16) = 1;
        localbd(3 * ny / 4 + 1 + nx * 2 + ny:ny + nx * 2 + ny, 16) = -1;
    end
end
end
