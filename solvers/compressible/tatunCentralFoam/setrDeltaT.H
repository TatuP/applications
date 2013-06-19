{       
	        
	    // Set the reciprocal time-step from the local Courant number
    	/*rDeltaT.dimensionedInternalField() = max
	(   
			1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
			 fvc::surfaceSum(amaxSf)().dimensionedInternalField()
			/(maxCo*mesh.V())
	); */ 
	surfaceScalarField amaxSfbyDelta
	(
	 mesh.surfaceInterpolation::deltaCoeffs()*amaxSf
	);

	rDeltaT.dimensionedInternalField() = max
	(
	 1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
	 fvc::surfaceSum(mag(amaxSfbyDelta/mesh.magSf()))/(2.0*maxCo)
	); 
	// Update tho boundary values of the reciprocal time-step
	rDeltaT.correctBoundaryConditions();
		
}
