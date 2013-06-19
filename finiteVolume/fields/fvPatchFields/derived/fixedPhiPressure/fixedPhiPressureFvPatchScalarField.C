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

#include "fixedPhiPressureFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedPhiPressureFvPatchScalarField::fixedPhiPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    UName_("U"),
    phiName_("phi"),
    rhoName_("rho"),
    adjoint_(false),
    flowRate_(0.0)
{}


Foam::fixedPhiPressureFvPatchScalarField::fixedPhiPressureFvPatchScalarField
(
    const fixedPhiPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    adjoint_(ptf.adjoint_),
    flowRate_(readScalar(dict.lookup("flowRate")))
    //flowRate_(ptf.flowRate_)
    //flowRate_(ptf.flowRate_().clone().ptr())
{}


Foam::fixedPhiPressureFvPatchScalarField::fixedPhiPressureFvPatchScalarField
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
    adjoint_(dict.lookup("adjoint")),
    flowRate_(DataEntry<scalar>::New("flowRate", dict))
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


Foam::fixedPhiPressureFvPatchScalarField::fixedPhiPressureFvPatchScalarField
(
    const fixedPhiPressureFvPatchScalarField& wbppsf
)
:
    fixedGradientFvPatchScalarField(wbppsf),
    UName_(wbppsf.UName_),
    phiName_(wbppsf.phiName_),
    rhoName_(wbppsf.rhoName_),
    adjoint_(wbppsf.adjoint_),
    flowRate_(wbppsf.flowRate_().clone().ptr())
{}


Foam::fixedPhiPressureFvPatchScalarField::fixedPhiPressureFvPatchScalarField
(
    const fixedPhiPressureFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    UName_(wbppsf.UName_),
    phiName_(wbppsf.phiName_),
    rhoName_(wbppsf.rhoName_),
    adjoint_(wbppsf.adjoint_),
    flowRate_(wbppsf.flowRate_().clone().ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedPhiPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& Up =
        patch().lookupPatchField<volVectorField, vector>(UName_);

    const surfaceScalarField& phi =
        db().lookupObject<surfaceScalarField>(phiName_);

    fvsPatchField<scalar> phip =
        patch().patchField<surfaceScalarField, scalar>(phi);
    
    //fvsPatchField<scalar> desiredPhi = flowRate_/patch().magSf(); 
    const fvPatchField<scalar>& rhop =
            patch().lookupPatchField<volScalarField, scalar>(rhoName_);
	
    if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
    {    

        phip /= rhop;
	
    }

    comst scalarField desiredU = flowRate_/rhop;

    const fvPatchField<scalar>& rAp =
        patch().lookupPatchField<volScalarField, scalar>("(1|A("+UName_+"))");


    const scalar t = db().time().timeOutputValue();
//   const fvPatchField<scalar>& rhop =
//        patch().lookupPatchField<volScalarField, scalar>(rhoName_);

    if (adjoint_)
    {
        gradient() = ((patch().Sf() & Up) - phip)/patch().magSf()/rAp;
    }
    else
    {
        //gradient() = ( - (patch().Sf() & Up))/patch().magSf()/rAp; ///rhop
    }
//    Info<< min(gradient()) << " < gradient() < " << max(gradient()) << endl;
//    Info<< min(rAp) << " < rAp < " << max(rAp) << endl;

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::fixedPhiPressureFvPatchScalarField::write(Ostream& os) const
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
        fixedPhiPressureFvPatchScalarField
    );
}

// ************************************************************************* //
