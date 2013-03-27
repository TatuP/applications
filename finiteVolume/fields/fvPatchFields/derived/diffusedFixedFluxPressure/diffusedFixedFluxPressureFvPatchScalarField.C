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

\*---------------------------------------------------------------------------*/

#include "diffusedFixedFluxPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diffusedFixedFluxPressureFvPatchScalarField::diffusedFixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    adjoint_(false)
{}


Foam::diffusedFixedFluxPressureFvPatchScalarField::diffusedFixedFluxPressureFvPatchScalarField
(
    const diffusedFixedFluxPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    adjoint_(ptf.adjoint_)
{}


Foam::diffusedFixedFluxPressureFvPatchScalarField::diffusedFixedFluxPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    UName_(dict.lookupOrDefault<word>("U", "U")),
    phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    rhoName_(dict.lookupOrDefault<word>("rho", "rho")),
    adjoint_(dict.lookup("adjoint"))
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
}


Foam::diffusedFixedFluxPressureFvPatchScalarField::diffusedFixedFluxPressureFvPatchScalarField
(
    const diffusedFixedFluxPressureFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    UName_(wbppsf.UName_),
    phiName_(wbppsf.phiName_),
    rhoName_(wbppsf.rhoName_),
    adjoint_(wbppsf.adjoint_)
{}


Foam::diffusedFixedFluxPressureFvPatchScalarField::diffusedFixedFluxPressureFvPatchScalarField
(
    const diffusedFixedFluxPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    UName_(wbppsf.UName_),
    phiName_(wbppsf.phiName_),
    rhoName_(wbppsf.rhoName_),
    adjoint_(wbppsf.adjoint_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diffusedFixedFluxPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar dp = 3.5;
    /*const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    fvsPatchField<scalar> phip =
        patch().patchField<surfaceScalarField, scalar>(phi);*/
  //  Info << "start of updateCoeffs" << endl;
    const surfaceScalarField& phiHbyA =
	        db().lookupObject<surfaceScalarField>("(1|A("+UName_+"_))");

    const surfaceScalarField& phi =
	    db().lookupObject<surfaceScalarField>(phiName_);

    fvsPatchField<scalar> phiHbyAp =
		patch().patchField<surfaceScalarField, scalar>(phiHbyA);

    fvsPatchField<scalar> phip =
	    patch().patchField<surfaceScalarField, scalar>(phi);


   /* if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {
        const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);

        phip /= rhop;
	
    }

    const fvPatchField<scalar>& rAp =
        patch().lookupPatchField<volScalarField, scalar>("(1|A("+UName_+"))");
*/
   /* const fvPatchField<scalar>& rAp =
        patch().lookupPatchField<volScalarField, scalar>("(1|A("+UName_+"))");*/

//   const fvPatchField<scalar>& rhop =
//        patch().lookupPatchField<volScalarField, scalar>(rhoName_);
//    Info << "before gradient() evaluation" << endl;
    if (adjoint_)
    {
        //gradient() = ((patch().Sf() & Up) - phip)/patch().magSf()/rAp/dp;
	    gradient() = (phip - phiHbyAp/patch().magSf())/patch().magSf();
    }
    else
    {
        //gradient() = (phip - (patch().Sf() & Up))/patch().magSf()/rAp/dp; ///rhop
	    gradient() = (phiHbyAp - phip)/patch().magSf()/10.0;
    }
//    Info<< min(gradient()) << " < gradient() < " << max(gradient()) << endl;
//    Info<< min(rAp) << " < rAp < " << max(rAp) << endl;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::diffusedFixedFluxPressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntryIfDifferent<word>(os, "U", "U", UName_);
    writeEntryIfDifferent<word>(os, "phi", "phi", phiName_);
    writeEntryIfDifferent<word>(os, "rho", "rho", rhoName_);
    os.writeKeyword("adjoint") << adjoint_ << token::END_STATEMENT << nl;
    gradient().writeEntry("gradient", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        diffusedFixedFluxPressureFvPatchScalarField
    );
}

// ************************************************************************* //
