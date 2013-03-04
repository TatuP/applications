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

    //thermo.correct();
    combustion->correct();
    dQ = combustion->dQ();
    label inertIndex = -1;
    volScalarField Yt = 0.0*Y[0];

    for (label i=0; i<Y.size(); i++)
    {
        if (Y[i].name() != inertSpecie)
        {
            volScalarField& Yi = Y[i];
	    volScalarField& rhoYi = rhoY[i]; 
//-------Inviscid predictor + viscous corrector for Yi begins-----------


	    // interpolation of the Yi to the cell faces
            surfaceScalarField Yi_pos =
    	        fvc::interpolate(Yi, pos, "reconstruct(Yi)"); // interpolate Yi at "positive" face
            surfaceScalarField Yi_neg =
                fvc::interpolate(Yi, neg, "reconstruct(Yi)"); // interpolate Yi at "negative" face

	    // weighted average of phi*Yi based on local speed of sound at the "positive" and "negative" cell faces,
	    // where aphiv_pos and aphiv_neg are the phi values at the cell faces, multiplied with the weighting factors 
	    // calculated based on local speed of sound at the cell faces
	    surfaceScalarField phiYi("phiYi", aphiv_pos*rho_pos*Yi_pos + aphiv_neg*rho_neg*Yi_neg); // phiYi=rho*U*Yi
	    
	    // ideally there should be a predictor value for the specie source, rho*combustion->R(Yi)
	    // however, I have not been able to modify the combustion->R function to support the solution of field rhoYi
	    /*volScalarField RYiPredi 
	    (
		    IOobject
		    (
			"RYiPredi",
			runTime.timeName(),
		        mesh,
        		IOobject::NO_READ,
		        IOobject::NO_WRITE
		    ),
	    rho*combustion->R(Yi)
	    );*/

	    
            solve
            (
		fvm::ddt(rhoYi)
	      + fvc::div(phiYi)
              //== combustion->Rrho(Yi,rhoYi)
            );
	    Yi = rhoYi/rho;
	    Yi.correctBoundaryConditions();
	    Yi.max(0.0);

            rhoYi.boundaryField() = rho.boundaryField()*Yi.boundaryField();

	    // Solve the viscous corrector for Yi (if the predictor is used)
            solve
            (
		fvm::ddt(rho,Yi) - fvc::ddt(rho,Yi)
	      + fvm::div(phi,Yi) - fvc::div(phi,Yi)
	      - fvm::laplacian(muEff/0.7,Yi) // turbulent Schmidt number 1.0 at the moment 
	        == combustion->R(Yi) // - RYiPredi/rho (in an ideal case)
	      ,
		mesh.solver("Yi")	
            ); 


// -----------Inviscid predictor + viscous corrector for Yi ends-----------

// ---------reactingFoam Yi eqn---------- (not used)
	    // species transport equation exactly as in reactingFoam,
	    // where it is assumed that Schmidt number Sc = 1.0, so that
	    // effective mu (turbulent + laminar mu) is exactly the diffusion
	    // coefficient in the specie transport equation
	/*    solve
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              - fvm::laplacian(turbulence->muEff(), Yi)
             == combustion->R(Yi)
	     ,
                mesh.solver("Yi")
            );*/


	    Yi.correctBoundaryConditions();
	    // print the mass fraction value range, if it is nonphysical (for debugging)
	    scalar minYi = min(Yi).value();
	    scalar maxYi = max(Yi).value();
            if ((minYi < -1e-100) || (maxYi > 1.0)) // print Yi range when it is unphysical
            {
    		Info<< minYi << " < Yi < " << maxYi << endl;
	    }
            Yi.max(0.0);
	    Yi.min(1.0);
	    rhoYi = rho*Yi; // store rhoYi for predictor equation of the next time step
	    
            Yt += Yi;
        }
        else
        {
            inertIndex = i;
        }
    }
 
    Y[inertIndex] = scalar(1) - Yt;
    volScalarField& Yin = Y[inertIndex];

// ---debugging stuff starts---
    // print the mass fraction value range of the inert specie, if it is unphysical
    scalar minYin = min(Yin).value();
    scalar maxYin = max(Yin).value();
    if ((minYin < 0.0) || (maxYin > 1.0))
    {
    	Info<< minYin << " < Yinert < " << maxYin << endl;
    }
    scalar minYsum = min(Yt + Yin).value();
    scalar maxYsum = max(Yt + Yin).value();
    // print the sum of each mass fractions if it is unphysical (should be exactly 1.0)
    if ((minYsum < 1.0) || (maxYsum > 1.0))
    {
	    Info << minYsum << " < sum(Yi) < " << maxYsum << endl;
    }
    YinertNeg = Yin; // YinertNeg is an extra field which shows negative mass fractions (for debugging)
    YinertNeg.min(0.0); // this extra field should have a maximum value of 0
// ---debugging stuff ends---
    
    Y[inertIndex].max(0.0);

    
}
