	scalar targetPhiF = 0.1*9.1980e-5;//9.198e-5; // target phi value per cell for fuelAir channel (two cells in y direction)
	scalar targetPhiCo = 0.1*5.2791e-5;//5.2791e-5;// target phi for cooler channel (upper)
	scalar targetPhiCa = 0.7e-5;//2.28e-6;// target phi for carrier channel (lower)
	scalar targetRho = 5.0;

	dimensionedScalar tau("tau",dimensionSet(0,0,1,0,0,0,0),8.0*runTime.deltaT().value()); // relaxation time scale
	scalar xRealInlet = -1.56e-2; // position of the "real" inlet, where the converging part of the nozzle starts
	scalar yMinF = 0.0008; // 0.0010; // y-range where fuelAir channel is 
	scalar yMaxF = 0.00282;// 

	volScalarField x = mesh.C().component(0);
	volScalarField y = mesh.C().component(1);
	scalar xMin = min(x).value(); // minimum value in the mesh (where the actual inlet is)
	scalar xMax = xMin + 0.8*(xRealInlet-xMin); // range where the phi correction is done (30 % from the inlet)
        //#include "targetRho.H"
	diffRhoU = rho*U*0.0; // initialize the difference between current rho*U and target rho*U (in createFields.H)
	volVectorField sigmaDiffRhoU = diffRhoU; // sigma*diffRhoU
	//volScalarField sigmaDiffRho(rho*0.0);
	forAll(x,i) // loop through all cells
	{
		scalar xi = x[i]; // x position of cell i
		scalar AcellFace = 0.02e-2*(0.2004695-0.0919226)*1e-2; //mesh.magSf()[i]; // cell face area (in fuelAir channel)
		rhoUA[i].component(0) = rho[i]*U[i].component(0)*AcellFace;
		if (xi < xMax) // consider cells that are on the left side of xMax
		{
			scalar a = 10.0*(xi-xMin)/(xMax-xMin); // cell position scaled to 0 <= a <= 6.0
			scalar sigma = 0.5*(Foam::tanh(-a+7.0)+1.0);  // 0 < a < 6: sigma(a=0) ~ 1, sigma(a=6) ~ 0

			scalar diffRhoUi = rho[i]*U[i].component(0); // substract the target rho*Ux from this value

			if ( (y[i] < yMaxF ) && (y[i] > yMinF) ) // restrict to cells in fuelAir channel (mid)
			{
				diffRhoUi -= targetPhiF/AcellFace;
			}
			else if (y[i] > yMaxF) // cells in cooler channel (upper)
			{
				diffRhoUi -= targetPhiCo/AcellFace;
			}
			else // cells in carrier channel (lower)
			{
				diffRhoUi -= targetPhiCa/AcellFace;
			}
			diffRhoU[i].component(0) = diffRhoUi; // difference between current and target rho*Ux
			sigmaDiffRhoU[i].component(0) = sigma*diffRhoUi; // the weighted difference
			//sigmaDiffRho[i] = sigma*(rho[i] - targetRho);
		}
	}
	//tauMC = muEff*dev2(Foam::T(fvc::grad(U)));
	solve
	(
		  fvm::ddt(rho,U) - fvc::ddt(rho,U) // substract values computed earlier in the solver
		+ fvm::div(phi,U) - fvc::div(phi,U) 
		//- fvm::laplacian(muEff,U) - fvc::ddt(muEff,U)
		//- fvc::div(tauMC)
		//+ fvc::grad(p)		
	     	== -sigmaDiffRhoU/tau
	);

	rhoU = rho*U;
