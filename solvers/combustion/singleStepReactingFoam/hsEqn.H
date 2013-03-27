{
    fvScalarMatrix hsEqn
    (
        fvm::ddt(rho, hs)
      + mvConvection->fvmDiv(phi, hs)
      - fvm::laplacian(alphaEff,hs)//- fvm::laplacian(turbulence->alphaEff(), hs) // alphaEff calculated in sarkarCorrection
//      - fvm::laplacian(turbulence->muEff(), hs)  // unit lewis no.
     ==
        dpdt
      - (fvc::ddt(rho, K) + fvc::div(phi, K))
      + combustion->Sh()
    );

    hsEqn.relax();
    hsEqn.solve();

    thermo.correct();

    Info<< "T gas min/max   = " << min(T).value() << ", "
        << max(T).value() << endl;
}
