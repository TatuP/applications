/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rhoCentralFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
//#include "basicPsiThermo.H"
#include "turbulenceModel.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedRhoFvPatchScalarField.H"

//#include "hsCombustionThermo.H" //nakul
//#include "psiChemistryModel.H" //nakul
//#include "hsPsiMixtureThermo.H"
//#include "psiChemistryCombustionModel.H" // Tatu 
//#include "rhoChemistryCombustionModel.H" // tatu
#include "psiCombustionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"
    #include "readThermophysicalProperties.H"
    #include "readTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // --- upwind interpolation of primitive fields on faces

        surfaceScalarField rho_pos
        (
            fvc::interpolate(rho, pos, "reconstruct(rho)")
        );
        surfaceScalarField rho_neg
        (
            fvc::interpolate(rho, neg, "reconstruct(rho)")
        );

        surfaceVectorField rhoU_pos
        (
            fvc::interpolate(rhoU, pos, "reconstruct(U)")
        );
        surfaceVectorField rhoU_neg
        (
            fvc::interpolate(rhoU, neg, "reconstruct(U)")
        );

        volScalarField rPsi(1.0/psi);
        surfaceScalarField rPsi_pos
        (
            fvc::interpolate(rPsi, pos, "reconstruct(T)")
        );
        surfaceScalarField rPsi_neg
        (
            fvc::interpolate(rPsi, neg, "reconstruct(T)")
        );
// tatu start: internal energy e --> sensible enthalpy hs = e + p
/*        surfaceScalarField e_pos
        (
            fvc::interpolate(e, pos, "reconstruct(T)")
        );
        surfaceScalarField e_neg
        (
            fvc::interpolate(e, neg, "reconstruct(T)")
        );*/

        surfaceScalarField hs_pos
        (
            fvc::interpolate(hs, pos, "reconstruct(T)")
        );
        surfaceScalarField hs_neg
        (
            fvc::interpolate(hs, neg, "reconstruct(T)")
        );
// tatu end
        surfaceVectorField U_pos(rhoU_pos/rho_pos);
        surfaceVectorField U_neg(rhoU_neg/rho_neg);

        surfaceScalarField p_pos(rho_pos*rPsi_pos);
        surfaceScalarField p_neg(rho_neg*rPsi_neg);

        surfaceScalarField phiv_pos(U_pos & mesh.Sf());
        surfaceScalarField phiv_neg(U_neg & mesh.Sf());

        //volScalarField c(sqrt(thermo.Cp()/(thermo.Cp()-R)*rPsi));
	volScalarField c(sqrt(thermo.Cp()/(thermo.Cp()-rPsi/T)*rPsi));
	//volScalarField c1(sqrt(thermo.Cp()/(thermo.Cv())*rPsi));
	//volScalarField c(sqrt(thermo.gamma()/rPsi));
        //Info << min(c-c1).value() << "c-c1" << max(c-c1).value() << endl;
	soundSpeed = c;
	Mach = mag(rhoU)/(c*rho);
        surfaceScalarField cSf_pos
        (
            fvc::interpolate(c, pos, "reconstruct(T)")*mesh.magSf()
        );
        surfaceScalarField cSf_neg
        (
            fvc::interpolate(c, neg, "reconstruct(T)")*mesh.magSf()
        );

        surfaceScalarField ap
        (
            max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero)
        );
        surfaceScalarField am
        (
            min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero)
        );

        surfaceScalarField a_pos(ap/(ap - am));

        surfaceScalarField amaxSf("amaxSf", max(mag(am), mag(ap)));

        surfaceScalarField aSf(am*a_pos);

        if (fluxScheme == "Tadmor")
        {
            aSf = -0.5*amaxSf;
            a_pos = 0.5;
        }

        surfaceScalarField a_neg(1.0 - a_pos);

        phiv_pos *= a_pos;
        phiv_neg *= a_neg;

        surfaceScalarField aphiv_pos(phiv_pos - aSf);
        surfaceScalarField aphiv_neg(phiv_neg + aSf);

        // Reuse amaxSf for the maximum positive and negative fluxes
        // estimated by the central scheme
        amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));

        #include "compressibleCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        phi = aphiv_pos*rho_pos + aphiv_neg*rho_neg;

        surfaceVectorField phiUp
        (
            (aphiv_pos*rhoU_pos + aphiv_neg*rhoU_neg)
          + (a_pos*p_pos + a_neg*p_neg)*mesh.Sf()
        );
// tatu start: e --> hs (change internal energy to sensible enthalpy, to handle reactions)
        /*surfaceScalarField phiEp
        (
            aphiv_pos*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
          + aphiv_neg*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
          + aSf*p_pos - aSf*p_neg
        );*/

	surfaceScalarField phiHp
        (
            aphiv_pos*(rho_pos*(hs_pos + 0.5*magSqr(U_pos))) //+ p_pos)
          + aphiv_neg*(rho_neg*(hs_neg + 0.5*magSqr(U_neg))) //+ p_neg)
          + aSf*p_pos - aSf*p_neg
        );
// tatu end
        volScalarField muEff(turbulence->muEff());
	volScalarField mut(turbulence->mut());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.dimensionedInternalField() =
            rhoU.dimensionedInternalField()
           /rho.dimensionedInternalField();
        U.correctBoundaryConditions();
        rhoU.boundaryField() = rho.boundaryField()*U.boundaryField();

        volScalarField rhoBydt(rho/runTime.deltaT());


	
        if (!inviscid)
        {
            //solve
	    tmp<fvVectorMatrix> UEqn
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
	    rUA = 1.0/UEqn().A();  // store A matrix for fixedFluxPressure BC
	    solve(UEqn());
            rhoU = rho*U;
        }
// tatu start: add specie transport equation, replace rhoE --> rhoH and e --> hs
	#include "YEqn.H" //nakul

        // --- Solve energy
        surfaceScalarField sigmaDotU
        (
            (
                fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
              + (mesh.Sf() & fvc::interpolate(tauMC))
            )
            & (a_pos*U_pos + a_neg*U_neg)
        );

	volScalarField dpdt = fvc::ddt(p);

	volScalarField shPredi = combustion->Sh();
        solve
        (
	    fvm::ddt(rhoH)
          + fvc::div(phiHp)
          - fvc::div(sigmaDotU)
	  - dpdt // rhoE = rhoH - p
	  == shPredi // include enthalpy calculation in the inviscid predictor equation?
        );

        //e = rhoE/rho - 0.5*magSqr(U);
	hs = rhoH/rho - 0.5*magSqr(U);
        hs.correctBoundaryConditions();//e.correctBoundaryConditions();
        thermo.correct();
        rhoH.boundaryField() = // rhoE
            rho.boundaryField()*
            (
                hs.boundaryField() + 0.5*magSqr(U.boundaryField())
            );

	//volScalarField shDiff = combustion->Sh() - shPredi;
        //Info << min(shDiff).value() << " < shDiff < " << max(shDiff).value() << endl;
	//volScalarField TrPsi = T+rPsi;
        if (!inviscid)
        {
            volScalarField k("k", thermo.Cp()*muEff/Pr);//thermo.Cp()*muEff/Pr);
            solve
            (
                fvm::ddt(rho, hs) - fvc::ddt(rho, hs)//fvm::ddt(rho, e) - fvc::ddt(rho, e)
              - fvm::laplacian(turbulence->alphaEff(), hs) // alphaEff = alpha + alphat 
	      + fvc::laplacian(turbulence->alphaEff(), hs)  // "remove" the contribution from the inviscid predictor
              - fvc::laplacian(k, T)   // originally minus sign!
	      //== combustion->Sh()// - shPredi
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
	
	//#include "/home/tptatu/OpenFOAM/tptatu-2.1.1/applications/targetFields/targetPhi.H"

        /*p.dimensionedInternalField() =
            rho.dimensionedInternalField()
           /psi.dimensionedInternalField();
        p.correctBoundaryConditions();
        rho.boundaryField() = psi.boundaryField()*p.boundaryField();*/

	/*rhoU = rho*U;
	rhoH = rho*(hs + 0.5*magSqr(U));*/
	K = 0.5*magSqr(U);
	h = thermo.Cp()*T + K;

        runTime.write();


	Info<< min(rho).value() <<" < rho < " << max(rho).value() << endl; // for debugging
	Info<< min(T).value() <<" < T < " << max(T).value() << endl;	   // for debugging
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
