	scalar targetPhiF = 9.1980e-5;//9.198e-5; // target phi value per cell for fuelAir channel (two cells in y direction)
	scalar targetPhiCo = 5.2791e-5;//5.2791e-5;// target phi for cooler channel (upper)
	scalar targetPhiCa = 0.7e-5;//2.28e-6;// Target phi for carrier channel (lower)
	scalar targetRho = 5.1;
	//dimensionedScalar targetT("targetT",dimensionSet(0,0,0,1,0,0,0),320.0);
	scalar targetT = 320.0;
	//dimensionedScalar unitT("unitT",dimensionSet(0,0,0,1,0,0,0),1.0);

	dimensionedScalar tau("tau",dimensionSet(0,0,1,0,0,0,0),8.0*runTime.deltaT().value()); // relaxation time scale
	scalar xRealInlet = -1.56e-2; // position of the "real" inlet, where the converging part of the nozzle starts
	scalar yMinF = 0.00067; // 0.0010; // y-range where fuelAir channel is 
	scalar yMaxF = 0.00242;// 

	volScalarField x = mesh.C().component(0);
	volScalarField y = mesh.C().component(1);
	scalar xMin = min(x).value(); // minimum value in the mesh (where the actual inlet is)
	scalar xMax = xMin + 0.9*(xRealInlet-xMin); // range where the phi correction is done (30 % from the inlet)
	xMin = xMin + 0.1*(xMax-xMin);
        //#include "targetRho.H"
	diffRhoU = rho*U*0.0; // initialize the difference between current rho*U and target rho*U (in createFields.H)
	volVectorField sigmaDiffRhoU = diffRhoU; // sigma*diffRhoU
	volScalarField sigmaDiffT = T*0.0;
	rhoUA = diffRhoU;
	//volScalarField sigmaDiffRho(rho*0.0);
	//Info << "before cell loop" << endl;
	//#include "targetRho.H"
	//thermo.correct();
	scalar AcellFace = 0.02e-2*(0.2004695-0.0919226)*1e-2; //mesh.magSf()[i]; // cell face area (in fuelAir channel)
	forAll(x,i) // loop through all cells
	{
		scalar xi = x[i]; // x position of cell i

		rhoUA[i].component(0) = rho[i]*U[i].component(0)*AcellFace;
		if ( (xi < xMax) && (xi > xMin) ) // consider cells that are on the left side of xMax
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
			diffT[i] = (T[i]-targetT);

			sigmaDiffT[i] = sigma*diffT[i];
			sigmaDiffRhoU[i].component(0) = sigma*diffRhoUi; // the weighted difference
			//sigmaDiffRho[i] = sigma*(rho[i] - targetRho);
		}
	}
	volScalarField Kold = 0.5*magSqr(U);
	solve
	(
		  fvm::ddt(rho,U) + fvm::div(phi,U) // substract values computed earlier in the solver
		- fvc::ddt(rho,U) - fvc::div(phi,U) 		
	     	== -sigmaDiffRhoU/tau
	);
	/*volScalarField Knew = 0.5*magSqr(U);
        solve
        (
                fvm::ddt(rho, hs) + fvc::ddt(rho,Knew) 
              - fvc::ddt(rho, hs) - fvc::ddt(rho,Kold)
	      == -rho*thermo.Cp()*sigmaDiffT/tau
        );
        thermo.correct();*/

        //rhoH = rho*(hs + 0.5*magSqr(U));
	//rhoU = rho*U;
