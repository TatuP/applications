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

Class
    Foam::massFlowTotalPressureFvPatchScalarField

Description
    Foam::massFlowTotalPressureFvPatchScalarField

SourceFiles
    massFlowTotalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef massFlowTotalPressureFvPatchScalarField_H
#define massFlowTotalPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class massFlowTotalPressureFvPatch Declaration
\*---------------------------------------------------------------------------*/

class massFlowTotalPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;

        //- Name of the compressibility field used to calculate the wave speed
        word psiName_;

        //- Heat capacity ratio
        scalar gamma_;

        //- Total pressure
        scalarField p0_;

	scalar massFlow_;//
	scalar Mave_;//



public:

    //- Runtime type information
    TypeName("massFlowTotalPressure");


    // Constructors

        //- Construct from patch and internal field
        massFlowTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        massFlowTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given massFlowTotalPressureFvPatchScalarField
        //  onto a new patch
        massFlowTotalPressureFvPatchScalarField
        (
            const massFlowTotalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        massFlowTotalPressureFvPatchScalarField
        (
            const massFlowTotalPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new massFlowTotalPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        massFlowTotalPressureFvPatchScalarField
        (
            const massFlowTotalPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new massFlowTotalPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& UName()
            {
                return UName_;
            }

            //- Return the name of the flux field
            const word& phiName() const
            {
                return phiName_;
            }

            //- Return reference to the name of the flux field
            //  to allow adjustment
            word& phiName()
            {
                return phiName_;
            }

            //- Return the name of the density field
            const word& rhoName() const
            {
                return rhoName_;
            }

            //- Return reference to the name of the density field
            //  to allow adjustment
            word& rhoName()
            {
                return rhoName_;
            }

            //- Return the name of the compressibility field
            const word& psiName() const
            {
                return psiName_;
            }

            //- Return reference to the name of the compressibility field
            //  to allow adjustment
            word& psiName()
            {
                return psiName_;
            }

            //- Return the heat capacity ratio
            scalar gamma() const
            {
                return gamma_;
            }

            //- Return reference to the heat capacity ratio to allow adjustment
            scalar& gamma()
            {
                return gamma_;
            }

            //- Return the total pressure
            const scalarField& p0() const
            {
                return p0_;
            }

            //- Return reference to the total pressure to allow adjustment
            scalarField& p0()
            {
                return p0_;
            }

            scalar massFlow() const
            {
                return massFlow_;
            }

            //- Return reference to the heat capacity ratio to allow adjustment
            scalar& massFlow()
            {
                return massFlow_;
            }

            scalar Mave() const
            {
                return Mave_;
            }

            //- Return reference to the heat capacity ratio to allow adjustment
            scalar& Mave()
            {
                return Mave_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            //  using the given patch total pressure and velocity fields
            virtual void updateCoeffs
            (
                const scalarField& p0p,
                const vectorField& Up
            );

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
