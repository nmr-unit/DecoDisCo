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

if strcmp(scan.type, 'ge-epi')
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
elseif strcmp(scan.type, 'se-epi')
    t2star = 1/(1/t2 + 1/t2prime);
    t2diff = 1/(1/t2 - 1/t2prime);
    if strcmp(scan.direction,'up')
        if strcmp(scan.pftype, 'zerofill')
            qplusdiff  = Ta/2/t2diff - 1i*pi*(B0*Ta + ycoord);
            qplusstar  = Ta/2/t2star - 1i*pi*(B0*Ta + ycoord);        
            PSF    = exp(1i*2*pi*B0*Te)*((exp((2*pf-1)*qplusdiff)-1)./(2*yres*qplusdiff) + (1-exp(-qplusstar)))./(2*yres*qplusstar);
        elseif strcmp(scan.pftype, 'conjfill')
            qplusdiff  = Ta/2/t2diff - 1i*pi*(B0*Ta + ycoord);
            qplusstar  = Ta/2/t2star - 1i*pi*(B0*Ta + ycoord);
            PSF = exp(1i*2*pi*B0*Te)*(((exp((2*pf-1)*qplusdiff)-1)./(2*yres*qplusdiff) + (1-exp(-qplusstar)))./(2*yres*qplusstar) + (exp((2*pf-1)*qplusdiff) - exp(qplusstar))./(2*yres*qplusstar));
        end
    elseif strcmp(scan.direction,'down')
        if strcmp(scan.pftype, 'zerofill')
            qminusdiff  = Ta/2/t2diff - 1i*pi*(B0*Ta - ycoord);
            qminusstar  = Ta/2/t2star - 1i*pi*(B0*Ta - ycoord); 
            PSF    = exp(1i*2*pi*B0*Te)*((exp((2*pf-1)*qminusdiff)-1)./(2*yres*qminusdiff) + (1-exp(-qminusstar)))./(2*yres*qminusstar);
        elseif strcmp(scan.pftype, 'conjfill')
            qminusdiff  = Ta/2/t2diff - 1i*pi*(B0*Ta - ycoord);
            qminusstar  = Ta/2/t2star - 1i*pi*(B0*Ta - ycoord);
            PSF = exp(1i*2*pi*B0*Te)*(((exp((2*pf-1)*qminusdiff)-1)./(2*yres*qminusdiff) + (1-exp(-qminusstar)))./(2*yres*qminusstar) + (exp((2*pf-1)*qminusdiff) - exp(qminusstar))./(2*yres*qminusstar));
        else
            error('Unknown Partial Fourier Type. Possible ''zerofill'' and ''conjfill''.');
        end
    else
        error('Unknown direction of traversal. Possible are ''up'' and ''down''.');
    end
elseif strcmp(scan.type, 'depicting')
    qplus  = -Ta/2/t2star + 1i*pi*(B0*Ta + ycoord);
    qminus = -Ta/2/t2star + 1i*pi*(B0*Ta - ycoord);
    PSF    = exp(1i*2*pi*Te*B0)*((1-exp(qplus))./qplus + (1-exp(qminus))./qminus)/2/yres;
else
    error('Unknown scan type. Possible are ''epi'' and ''depicting''.');
end

norm_factor = sqrt(trapz(real(PSF))^2 + trapz(imag(PSF))^2);
PSF = (real(PSF)/norm_factor + 1i*imag(PSF)/norm_factor);

end