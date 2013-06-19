/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License

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
    Foam::combustionModels::EDM

Description
    Simple infinitely fast chemistry combustion model based on the principle
    mixed is burnt. Additional parameter C is used to distribute the heat
    release rate.in time

SourceFiles
    EDM.C

\*---------------------------------------------------------------------------*/

#ifndef EDM_H
#define EDM_H
#include "singleStepCombustion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{

/*---------------------------------------------------------------------------*\
                  Class EDM Declaration
\*---------------------------------------------------------------------------*/

template<class CombThermoType, class ThermoType>
class EDM
:
    public singleStepCombustion <CombThermoType, ThermoType>
{
    // Private data

        //- Model constant
	scalar xIgnition_;
        scalar C_;
	Time& runTime_;
	//volScalarField rDeltaT_;


    // Private Member Functions

        //- Disallow copy construct
        EDM(const EDM&);

        //- Disallow default bitwise assignment
        void operator=(const EDM&);

//        singleStepReactingMixture<ThermoType>& singleMixture_;


public:

    //- Runtime type information
    TypeName("EDM");


    // Constructors

        //- Construct from components
        EDM
        (
            const word& modelType, const fvMesh& mesh
        );


    //- Destructor
    virtual ~EDM();


    // Member Functions

        // Evolution

            //- Correct combustion rate
            virtual void correct();


    // I-O

            //- Update properties
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "EDM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
