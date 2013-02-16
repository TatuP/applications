/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    tatunCentralSpeciesFoam

Description
    Viscous, compressible, reactive flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

 #include "fvCFD.H"
 #include "hsCombustionThermo.H" //nakul
 #include "psiChemistryModel.H" //nakul
 #include "hsPsiMixtureThermo.H"
 #include "turbulenceModel.H"    //Arv 
 #include "zeroGradientFvPatchFields.H"
 #include "fixedRhoFvPatchScalarField.H"
 #include "multivariateScheme.H" //nakul

 #include "psiChemistryCombustionModel.H" // Tatu 
 #include "rhoChemistryCombustionModel.H" // tatu
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readChemistryProperties.H"
    #include "createFields.H"
    #include "readThermophysicalProperties.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *//
    #include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero",dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // --- upwind interpolation of primitive fields on faces

        surfaceScalarField rho_pos =
            fvc::interpolate(rho, pos, "reconstruct(rho)");
        surfaceScalarField rho_neg =
            fvc::interpolate(rho, neg, "reconstruct(rho)");

        surfaceVectorField rhoU_pos =
            fvc::interpolate(rhoU, pos, "reconstruct(U)");
        surfaceVectorField rhoU_neg =
            fvc::interpolate(rhoU, neg, "reconstruct(U)");

        volScalarField rPsi = 1.0/psi; // reverse psi
        surfaceScalarField rPsi_pos =
            fvc::interpolate(rPsi, pos, "reconstruct(T)");
        surfaceScalarField rPsi_neg =
            fvc::interpolate(rPsi, neg, "reconstruct(T)");

        surfaceScalarField hs_pos =
            fvc::interpolate(hs, pos, "reconstruct(T)");
        surfaceScalarField hs_neg =
            fvc::interpolate(hs, neg, "reconstruct(T)");

        surfaceVectorField U_pos = rhoU_pos/rho_pos;
        surfaceVectorField U_neg = rhoU_neg/rho_neg;

        surfaceScalarField p_pos = rho_pos*rPsi_pos;
        surfaceScalarField p_neg = rho_neg*rPsi_neg;

        surfaceScalarField phiv_pos = U_pos & mesh.Sf(); // phi_f terms in eq. (8), unit [kg/m^2/s]
        surfaceScalarField phiv_neg = U_neg & mesh.Sf();

        volScalarField c = sqrt(thermo.Cp()/(thermo.Cp()-R)*rPsi); // local speed of sound, R is molar gas constant
//	volScalarField c = sqrt(thermo.Cp()/thermo.Cv()*rPsi); // local speed of sound, R is molar gas constant
        surfaceScalarField cSf_pos = fvc::interpolate(c, pos, "reconstruct(T)")*mesh.magSf(); // c_f*|S_f| term for (8)
        surfaceScalarField cSf_neg = fvc::interpolate(c, neg, "reconstruct(T)")*mesh.magSf();

        surfaceScalarField ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero); // volumetric flux coefficient associated with local 
        surfaceScalarField am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero); // speed of propagation. finishes eq. (8) 

        surfaceScalarField a_pos = ap/(ap - am); // Upwinding coefficient, eq. (9), based on local speed of sound (KNP method), acting as the limiter.
						 // this is gives +face contribution to a field value on the face

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap))); // w_f for upwinded method, eq. (10) in the article

        #include "compressibleCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        surfaceScalarField aSf = am*a_pos; // alpha*(1-alpha) for KNP method, used for diffusive flux upwinding coefficient w_f

        if (fluxScheme == "Tadmor") // no upwinding, symmetrical interpolation
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg = (1.0 - a_pos); // the -face contribution to a field value on the face, see eq. (9)

        phiv_pos *= a_pos; // phiv_pos, that is phi_f+ in article, multiplied with alpha, the upwinding coefficient
        phiv_neg *= a_neg; // phi_f-multiplied with (1-alpha)

        surfaceScalarField aphiv_pos = phiv_pos - aSf; // coefficient for a tensor field Psi_f+ and Psi_f- in reordered eq. (7).
        surfaceScalarField aphiv_neg = phiv_neg + aSf; // upwinded average phi_f*Psi_f = aphiv_pos*Psi_f+ + aphiv_neg*Psi_f-

        surfaceScalarField phi("phi", aphiv_pos*rho_pos + aphiv_neg*rho_neg); //phi=rho*U ; aphiv=U

        surfaceVectorField phiUp = // div(ur*ho_U) + nabla(p) for eq. (18), which is the inviscid contribution to the momentum equation
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg) // phiUp=rho*U*U + p
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf(); 

        surfaceScalarField phiHp =
            aphiv_pos*(rho_pos*(hs_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(hs_neg + 0.5*magSqr(U_neg)) + p_neg) //phiHp= U*(rho*Hstot + p)
          + aSf*p_pos - aSf*p_neg;

	volScalarField muEff(turbulence->muEff());
	volScalarField alphaEff(turbulence->alphaEff());

        volTensorField tauMC("tauMC", muEff*dev2(fvc::grad(U)().T())); // Deviatoric Part of Stress Tensor

        //#include "chemistry.H"
	
        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi)); // phi is scalar, defined on each cell faces, from which div is calculated using gauss's theorem
        Info<<"rho solved diagonally" << endl;
        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp)); // eq. (18), which is the inviscid contribution to the momentum equation (momentum predictor)
	Info<<"rhoU solved diagonally" << endl;

        U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField(); // This is done so that we can use conventional boundary conditions applied for U
        U.correctBoundaryConditions(); // interpolate U to boundary based on internal values?
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

        volScalarField rhoBydt(rho/runTime.deltaT()); // New variable rho/t

	

        if (!inviscid) // viscous contribution to the momentum flux, rhoU, from eq. (19)
        {
            solve // U solved (viscous momentum correction)
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
                - fvm::laplacian(muEff, U) //This formulation of Momentum Equation is kept unchanged for cold flow
              - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }
        #include "YEqn.H" //nakul
        
	// --- Solve energy
        surfaceScalarField sigmaDotU =
        (
            (
                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
              + (mesh.Sf() & fvc::interpolate(tauMC))
            )
            & (a_pos*U_pos + a_neg*U_neg)
        ); // (mu*grad(U) + Deviatoric Stress Tensor).U

        solve // inviscid energy predictor
        (
            fvm::ddt(rhoH)
          + fvc::div(phiHp)
          - fvc::div(sigmaDotU)
	  - fvc::ddt(p) //nakul
         == chemistrySh
	);
	Info<<"rhoH solved diagonally" << endl;
        hs = rhoH/rho - 0.5*magSqr(U);
        hs.correctBoundaryConditions();
        thermo.correct();
        rhoH.boundaryField() =
            rho.boundaryField()*
            (
                hs.boundaryField() + 0.5*magSqr(U.boundaryField())
            );

        if (!inviscid)
        {
            volScalarField k("k", thermo.Cp()*(turbulence->muEff()/Pr));
	    
            solve // diffusion equation for temperature corrector
            (
                fvm::ddt(rho, hs) - fvc::ddt(rho, hs)
              - fvm::laplacian(alphaEff, hs)
              + fvc::laplacian(alphaEff, hs)
              - fvc::laplacian(k, T) 
	      - fvc::ddt(p) //nakul
	     == chemistrySh
            );
            thermo.correct();
            rhoH = rho*(hs + 0.5*magSqr(U));
        }

        p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();
        p.correctBoundaryConditions();
        rho.boundaryField() = psi.boundaryField()*p.boundaryField();

         turbulence->correct();
         
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}
// ************************************************************************* //
