tmp<fv::convectionScheme<scalar> > mvConvection
(
    fv::convectionScheme<scalar>::New
    (
        mesh,
        fields,
        phi,
        mesh.divScheme("div(phi,Yi_h)")
    )
);

{
    combustion->correct();
    dQ = combustion->dQ();
    label inertIndex = -1;
    volScalarField Yt(0.0*Y[0]);

    forAll(Y, i)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff()/0.85, Yi)
             ==
                combustion->R(Yi)
            );

            YiEqn.relax();
            YiEqn.solve(mesh.solver("Yi"));

	    scalar YiMax = max(Yi).value();
	    scalar YiMin = min(Yi).value();
	    if ((YiMin < -1e-100) || (YiMax > 1.0))
	    {
		    Info << YiMin << " < Yi < " << YiMax << endl;
	    }
            Yi.max(0.0);
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }

    Y[inertIndex] = scalar(1) - Yt;
    volScalarField& Yin = Y[inertIndex];
    scalar YinMax = max(Yin).value();
    scalar YinMin = min(Yin).value();
    if ((YinMin < -1e-100) || (YinMax > 1.0))
    {
	    Info << YinMin << " < Yinert < " << YinMax << endl;
	    YinertNeg = Yin;
    }
    
    scalar YsumMax = max(Yt+Yin).value();
    scalar YsumMin = min(Yt+Yin).value();
    if ((YsumMin < 1.0) || (YsumMax > 1.0))
    {
	    Info<< YsumMin << " < sum(Yi) < " << YsumMax << endl;
    }
    Y[inertIndex].max(0.0);
}
