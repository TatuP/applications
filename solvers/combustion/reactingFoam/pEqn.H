rho = thermo.rho();

volScalarField rAU(1.0/UEqn.A());
U = rAU*UEqn.H();

if (pimple.transonic())
{
    surfaceScalarField phid
    (
        "phid",
        fvc::interpolate(psi)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        )
    );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvm::div(phid, p)
          - fvm::laplacian(rho*rAU, p)
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi == pEqn.flux();
        }
    }
}
else
{
    phi =
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf())
          + fvc::ddtPhiCorr(rAU, rho, U, phi)
        );

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::ddt(psi, p)
          + fvc::div(phi)
          - fvm::laplacian(rho*rAU, p)
        );

        pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pEqn.flux();
        }
    }
}

#include "rhoEqn.H"
#include "compressibleContinuityErrs.H"

Info<< "rho min/max : " << min(rho).value()
<< " " << max(rho).value() << endl;
Info<< "p min/max : " << min(p).value()
<< " " << max(p).value() << endl;
U -= rAU*fvc::grad(p);
U.correctBoundaryConditions();
K = 0.5*magSqr(U);

phiHbyA =
        fvc::interpolate(rho)
       *(
            (fvc::interpolate(U) & mesh.Sf()) //*(fvc::interpolate(asd) & mesh.Sf())
        );


dpdt = fvc::ddt(p);
