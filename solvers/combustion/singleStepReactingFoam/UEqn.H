 //fvVectorMatrix sarkarStressDiag =
 //     - fvm::laplacian(mut, U);
//volVectorField sarkarStressDev =
//      - fvc::div(mut*dev2(Foam::T(fvc::grad(U))));   
    fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(phi, U)
      //+ turbulence->divDevRhoReff(U)
      + sarkarStressDiag + sarkarStressDev
     ==
        rho*g
    );

    UEqn.relax();

    if (pimple.momentumPredictor())
    {
        solve(UEqn == -fvc::grad(p));
        K = 0.5*magSqr(U);
    }
