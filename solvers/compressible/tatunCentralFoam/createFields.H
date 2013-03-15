Info<< "Reading thermophysical properties\n" << endl;

// tatu begin: comment basicPsiThermo, and use either psi- or rho-based combustion chemsitry thermo

/*autoPtr<basicPsiThermo> pThermo
(
    basicPsiThermo::New(mesh)
);
basicPsiThermo& thermo = pThermo();*/


// rho-based
/*
autoPtr<combustionModels::rhoChemistryCombustionModel> combustion
(
    combustionModels::rhoChemistryCombustionModel::New
    (
        mesh
    )
); 
rhoChemistryModel& chemistry = combustion->pChemistry(); // rho-bsed
hsReactionThermo& thermo = chemistry.thermo(); // tatu added, for rhoThermo*/



/*// psi-based
autoPtr<combustionModels::psiChemistryCombustionModel> combustion
(
    combustionModels::psiChemistryCombustionModel::New
    (
        mesh
    )
); // psi-based
psiChemistryModel& chemistry = combustion->pChemistry(); // psi-based
hsCombustionThermo& thermo = chemistry.thermo(); */// tatu commented, for psiThermo

//------------------------------------------------------------------//

// psiCombustion model for singleStepCombustion (single step chemistry)
autoPtr<combustionModels::psiCombustionModel> combustion
(
    combustionModels::psiCombustionModel::New
    (
        mesh
    )
);
hsCombustionThermo& thermo = combustion->thermo();

basicMultiComponentMixture& composition = thermo.composition();

PtrList<volScalarField>& Y = composition.Y();

 
word inertSpecie(thermo.lookup("inertSpecie"));
// tatu done

volScalarField& p = thermo.p();
volScalarField& hs = thermo.hs(); // replaced e --> hs
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

/*forAll(Y,i)
{
	rhoY[i] = rho*Y[i];//dimensionedScalar("rhoYi", dimensionSet(1,-3,0,0,0), 1.0); //rho*Y[i];
	//rhoY[i]=dimensionSet(1,-3,0,0,0);
}*/

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    rho*U
);
// tatu start: rhoE --> rhoH
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

#include "compressibleCreatePhi.H"

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
combustion->setTurbulence(turbulence());
// tatu start: molar gas constant R because thermo.Cv is not available. needed to compute local speed of sound.
//		also initialize Y, hs etc. 		
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

PtrList<volScalarField> rhoY(Y.size());
forAll (Y, i)
{
         rhoY.set
         (
             i,
             new volScalarField
             (
                 IOobject
                 (
                     "rhoYi",
                     runTime.timeName(),
		     mesh,
                     IOobject::NO_READ,
                     IOobject::NO_WRITE
                 ),
                 rho*Y[i]
             )
         );
    	 fields.add(Y[i]);
}
fields.add(hs);

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
    //dimensionedScalar("zero", dimless, 0.0)
    dimensionedScalar("kappa", dimless, 0.0)
);
volScalarField rUA // rUA stored for fixedFluxPressure BC
(
 	 IOobject
 	 (
 	 "(1|A(U))",
	  runTime.timeName(),
	  mesh
	 ),
	 mesh,
	 //dimensionedScalar("rUA", dimTime, 1),
	 dimensionedScalar("rUA",dimensionSet(-1,3,1,0,0),1),
	 zeroGradientFvPatchScalarField::typeName
);
volScalarField dQ
(
    IOobject
    (
        "dQ",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("dQ", dimEnergy/dimTime, 0.0)
);
volScalarField YinertNeg // negative inert mass fraction field for debugging
(
    IOobject
    (
        "YinertNeg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("YinertNeg", dimEnergy/dimEnergy, 0.0)
);

volScalarField Mach
(
    IOobject
    (
        "Mach",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Mach", dimless, 0.0)
);
volScalarField soundSpeed
(
    IOobject
    (
        "soundSpeed",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mag(U)*0.0
);
volScalarField h 
(
    IOobject
    (
        "h",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0.5*magSqr(U)*0.0
);
volScalarField K 
(
    IOobject
    (
        "K",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0.5*magSqr(U)
);

volScalarField rDeltaT
(
    IOobject
    (
        "rDeltaT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    1/dimensionedScalar("maxDeltaT", dimTime, GREAT),
    zeroGradientFvPatchScalarField::typeName
);
// tatu end
