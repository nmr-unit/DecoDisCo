function PSF = getPSF(param, B0, t2star, scan)

Ny = param.ny;
yres = param.yres;
Te = param.Te;
Ta = param.Ta;
if isfield(scan, 'pf')
    pf = scan.pf;
end

PSF = zeros(Ny,1);
if mod(Ny,2) == 0
    ycoord = -Ny/2+1:Ny/2;
elseif mod(Ny,2) == 1
    ycoord   = -Ny/2+.5:Ny/2-.5;
end

if strcmp(scan.type, 'epi')
    if strcmp(scan.direction,'up')
        if strcmp(scan.pftype, 'zerofill')
            qplus  = -Ta/2/t2star + 1i*pi*(B0*Ta + ycoord);
            PSF    = exp(1i*2*pi*B0*Te)*(exp(-(2*pf-1)*qplus) - exp(qplus))./(2*yres*qplus);
        elseif strcmp(scan.pftype, 'conjfill')
            qplus  = -Ta/2/t2star + 1i*pi*(B0*Ta + ycoord);
            PSF = exp(1i*2*pi*B0*Te)*((exp(-(2*pf-1)*qplus) - exp(qplus))./(2*yres*qplus) + (exp((2*pf-1)*qplus) - exp(qplus))./(2*yres*qplus));
        end
    elseif strcmp(scan.direction,'down')
        if strcmp(scan.pftype, 'zerofill')
            qminus = -Ta/2/t2star + 1i*pi*(+B0*Ta - ycoord);
            PSF = exp(1i*2*pi*B0*Te)*(exp(-(2*pf-1)*qminus) - exp(qminus))./(2*yres*qminus);
        elseif strcmp(scan.pftype, 'conjfill')
            qminus = -Ta/2/t2star + 1i*pi*(+B0*Ta - ycoord);
            PSF = exp(1i*2*pi*B0*Te)*((exp(-(2*pf-1)*qminus) - exp(qminus))./(2*yres*qminus) + (exp((2*pf-1)*qminus) - exp(qminus))./(2*yres*qminus));
        else
            error('Unknown Partial Fourier Type. Possible ''zerofill'' and ''conjfill''.');
        end
    else
        error('Unknown direction of traversal. Possible are ''up'' and ''down''.');
    end
elseif strcmp(scan.type, 'depicting')
    qplus  = -Ta/2/t2star + 1i*pi*(B0*Ta + ycoord);
    qminus = -Ta/2/t2star + 1i*pi*(B0*Ta - ycoord);
    PSF(ic1) = exp(1i*2*pi*Te*B0)*((1-exp(qplus))./qplus + (1-exp(qminus))./qminus)/2/yres;
else
    error('Unknown scan type. Possible are ''epi'' and ''depicting''.');
end

norm_factor = sqrt(trapz(real(PSF))^2 + trapz(imag(PSF))^2);
PSF = (real(PSF)/norm_factor + 1i*imag(PSF)/norm_factor);

end