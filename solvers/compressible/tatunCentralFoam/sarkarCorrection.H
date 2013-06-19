volScalarField mut = turbulence->mut();
volScalarField alphat(turbulence->alphat());
volScalarField k(turbulence->k());

volScalarField epsilon(turbulence->epsilon());

//tmp<volTensorField> tgradU = fvc::grad(U);
//volScalarField G("RASModel::G", mut*(tgradU() && dev(twoSymm(tgradU()))));
//tgradU.clear();

tmp<volTensorField> tgradU = fvc::grad(U);
volScalarField S2((tgradU() && dev(twoSymm(tgradU()))));
tgradU.clear();


scalar eta0 = 4.38;
scalar beta = 0.012;
scalar Cmu = 0.0845;
volScalarField G("RASModel::G", mut*S2);

volScalarField eta(sqrt(mag(S2))*k/epsilon);
volScalarField eta3(eta*sqr(eta));
volScalarField R
(
	((eta*(-eta/eta0 + scalar(1)))/(beta*eta3 + scalar(1)))
);
epsilon.boundaryField().updateCoeffs();
volScalarField YM = 2.0*rho*epsilon*k/(soundSpeed*soundSpeed);
Info << min(YM).value() << " < YM < " << max(YM).value() << endl;

volScalarField kOld = k;
tmp<fvScalarMatrix> kEqn
(
        fvm::ddt(rho, k) //+ fvm::div(phi, k)
      - fvc::ddt(rho, k) //- fvc::div(phi, k)
	== -YM    
);


kEqn().relax();
solve(kEqn);
k.max(0.0);//bound(k,0);
Info << min(k-kOld).value() << " < k-kOld < " << max(k-kOld).value() << endl;

scalar C2 = 1.68;//1.92;
scalar C1 = 1.42;//1.44;


tmp<fvScalarMatrix> epsEqn
(
        fvm::ddt(rho, epsilon) //+ fvm::div(phi, epsilon)
      - fvc::ddt(rho, epsilon) //- fvc::div(phi, epsilon)
     ==
        (C1-R)*G*epsilon/k
      - (C1-R)*G*epsilon/kOld
      - fvm::Sp(C2*rho*epsilon/k, epsilon)
      + fvc::Sp(C2*rho*epsilon/k, epsilon)
);
epsEqn().relax();
epsEqn().boundaryManipulate(epsilon.boundaryField());

solve(epsEqn);
epsilon.max(0.0);

    // Re-calculate viscosity
volScalarField mutOld = mut;
mut = rho*Cmu*sqr(k)/epsilon;
Info << min(mut-mutOld).value() << " < mut-mutOld < " << max(mut-mutOld).value() << endl;
mut.correctBoundaryConditions();
    // Re-calculate thermal diffusivity
scalar Prt = 0.85;
alphat = mut/Prt;//mut_/Prt_;
alphat.correctBoundaryConditions();
volScalarField alphaEff = (turbulence->mu() + mut)/Prt;

fvVectorMatrix sarkarStressDiag =
      - fvm::laplacian(mut, U);
volVectorField sarkarStressDev =
      - fvc::div(mut*dev2(Foam::T(fvc::grad(U))));
