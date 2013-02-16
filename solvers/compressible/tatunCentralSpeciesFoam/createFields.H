Info<< "Reading thermophysical properties\n" << endl;

/*autoPtr<basicPsiThermo> pThermo
(
    basicPsiThermo::New(mesh)
);
basicPsiThermo& thermo = pThermo();*/

// tatu begin
autoPtr<combustionModels::rhoChemistryCombustionModel> combustion
(
    combustionModels::rhoChemistryCombustionModel::New
    (
        mesh
    )
); // rho-based
rhoChemistryModel& chemistry = combustion->pChemistry(); // rho-bsed
hsReactionThermo& thermo = chemistry.thermo(); // tatu added, for rhoThermo

/*
autoPtr<combustionModels::psiChemistryCombustionModel> combustion
(
    combustionModels::psiChemistryCombustionModel::New
    (
        mesh
    )
); // psi-based
psiChemistryModel& chemistry = combustion->pChemistry(); // psi-based
hsCombustionThermo& thermo = chemistry.thermo(); // tatu commented, for psiThermo
*/

basicMultiComponentMixture& composition = thermo.composition();

PtrList<volScalarField>& Y = composition.Y();
word inertSpecie(thermo.lookup("inertSpecie"));
// tatu done

volScalarField& p = thermo.p();
volScalarField& hs = thermo.hs();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();
const volScalarField& mu = thermo.mu();

bool inviscid(true);
if (max(mu.internalField()) > 0.0)
{
    inviscid = false;
}

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "rhoBoundaryTypes.H"
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho(),
    rhoBoundaryTypes
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);
// tatu start
volScalarField rhoH 
(
    IOobject
    (
        "rhoH",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(hs + 0.5*magSqr(U))
);
// tatu end
surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);


surfaceScalarField phi("phi", mesh.Sf() & fvc::interpolate(rhoU));

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
// tatu start
const volScalarField R
(
    IOobject
    (
        "R",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("R", dimensionSet(0,2,-2,-1,0), 8.3144)
);

multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll (Y, i)
{
    fields.add(Y[i]);
}
fields.add(hs);

//PtrList<volScalarField>& rhoY = rho*composition.Y();

/*forAll (Y, i)
{
    	rhoYi = rho*Y[i];

	fields.add(rhoYi);
}*/


DimensionedField<scalar, volMesh> chemistrySh
(
    IOobject
    (
        "chemistry::Sh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("chemistrySh", dimEnergy/dimTime/dimVolume, 0.0)
);

volScalarField kappa
(
    IOobject
    (
        "kappa",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);
// tatu end
