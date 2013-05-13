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

\*---------------------------------------------------------------------------*/

#include "EDM.H"
#include "wallFvPatch.H"
namespace Foam
{
namespace combustionModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
EDM<CombThermoType, ThermoType>::EDM
(
    const word& modelType, const fvMesh& mesh
)
:
    singleStepCombustion<CombThermoType, ThermoType>(modelType, mesh),
    C_(readScalar(this->coeffs().lookup("C")))
/*    singleMixture_
    (
        dynamic_cast<singleStepReactingMixture<ThermoType>&>(this->thermo())
    )*/
{}


// * * * * * * * * * * * * * * * * Destructors * * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
EDM<CombThermoType, ThermoType>::~EDM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CombThermoType, class ThermoType>
void EDM<CombThermoType, ThermoType>::correct()
{
    this->wFuel_ ==
        dimensionedScalar("zero", dimMass/pow3(dimLength)/dimTime, 0.0);
// tatu start

    tmp<volScalarField> trho(this->rho());
    const volScalarField& rho = trho();
    tmp<volScalarField> tepsilon(this->turbulence().epsilon());
    const volScalarField& epsilon = tepsilon();
    tmp<volScalarField> tk(this->turbulence().k());
    const volScalarField& k = tk();
    tmp<volScalarField> tmuEff(this->turbulence().muEff());
    const volScalarField& muEff = tmuEff();

    tmp<volScalarField> tx(this->mesh().C().component(0));
    const volScalarField x = tx();
    //const volScalarField& invTauTurb = epsilon/k;

    //const volScalarField& tau = this->turbulence().epsilon()/this->turbulence().k();
/*    const volScalarField& muEff = this->turbulence().muEff();
    const volScalarField& epsilon = this->turbulence().epsilon();
    const volScalarField& k = this->turbulence().k();*/

//    volScalarField& t = invTauTurb.max(0.001);

/*    scalar minTau = 1.0e7;
    scalar taui = minTau;
    forAll(tau,i)
    {
	
	taui = tau[i];
	if(taui < minTau) {
		minTau = taui;
        }
    }*/
    
//    Info<< "min(tau) = " << minTau << endl;
    scalar A = 4.0;
    scalar B = 0.5;

    if (this->active())
    {
       // this->singleMixture_.fresCorrect(); 

        const label fuelI = this->singleMixture_.fuelIndex();
        //const label O2I = this->singleMixture_.O2Inde
        const volScalarField& YFuel = this->thermo_->composition().Y()[fuelI];

        const dimensionedScalar s = this->singleMixture_.s();

        if (this->thermo_->composition().contains("O2"))
        {

	    //const Reaction<ThermoType>& reaction = this->operator[](0);	
            const volScalarField& YO2 = this->thermo_->composition().Y("O2");
	    //volScalarField YPsum = 0.0*YO2;
	    //this->singleMixture_.calcYprod();
	    //YPsum = this->singl
	    
	    //const List<scalar>& stoichList = this->singleMixture_.specieStoichCoeffs(); 
	    /*scalar stoichMPsum = 0.0; // sum(stoich coeff * molar mass) 
	    forAll(this->thermo_->composition().Y(),Yi)
	    {
                YPsum += Yi;	
		stoichMPsum += 0.0;*/
		//label indexi = Yi.index();
		/*if ((Yi.index() != YO2.index()) && (Yi != YFuel))
		{
			YPsum += Yi;
		}
            }*/
	    /*YPsum -= YFuel;
	    YPsum -= YO2;*/
	    /*Info<< min(YO2).value() << " < YO2 < " << max(YO2).value() << endl;
	    Info<< min(YFuel).value() << " < YFuel < " << max(YFuel).value() << endl;
	    Info<< min(epsilon).value() << " < epsilon < " << max(epsilon).value() << endl;
	    Info<< min(k).value() << " < k < " << max(k).value() << endl;
	    Info<< min(muEff).value() << " < muEff < " << max(muEff).value() << endl;*/

            /*this->wFuel_ ==
                 this->rho()/(this->mesh().time().deltaT()*C_)
                *min(YFuel, YO2/s.value());*/
	    //const volScalarField& rho = this->rho();
	    //Info<< min(rho).value() << " < rho < " << max(rho).value() << endl;
//	    scalar& wFuels = 0.0*rho;
	    scalar maxAtauimu = -1.0;
	    scalar minAtauimu = 1e15;
	    scalar maxAtauik = -1.0;
	    scalar minAtauik = 1e15;
	        //Info<< "asd" << endl;

	    forAll(rho,i) 
	    {
		scalar rhoi = rho[i];
		scalar Atauik = A*epsilon[i]/(k[i]+SMALL);
		scalar muEffi = muEff[i];
		scalar epsiloni = epsilon[i];
	        scalar tk =
                        0.1
                       *Foam::sqrt(muEffi/(rhoi+SMALL)/(epsiloni + SMALL));
		//maxtk = max(maxtk,tk);
		//mintk = min(mintk,tk);
	        scalar Atauimu = A/(tk+SMALL);
		maxAtauimu = max(maxAtauimu,Atauimu);
		minAtauimu = min(minAtauimu,Atauimu);
		maxAtauik = max(maxAtauik,Atauik);
		minAtauik = min(minAtauik,Atauik);
		//scalar Ataui = min(Atauik,Atauimu);
		scalar Ataui = Atauimu;
		scalar si = s.value();
		scalar YFueli = YFuel[i];
		scalar YO2si = YO2[i]/si; // s.value = stoichiometric fuel-oxygen ratio	
		scalar YProdi = 1.0-YFueli-YO2[i];
		scalar dt = this->mesh().time().deltaT().value();
		//Info<< "s.value() = " << s.value() << endl;
		/*if (Ataui>1.0/(dt*C_)) {
			Info<< "Ataui = " << Ataui << " > 1.0/(dt*C_) = " << 1.0/(dt*C_) << endl;
		}*/
		scalar YminReactants = min(YFueli,YO2si);
            	this->wFuel_[i] =
			min(Ataui,1.0/(dt*C_))*rho[i]*min(YminReactants,B*YProdi/(1.0+si));
			//min(Ataui,1.0/(dt*C_))*rho[i]*min(YminReactants, B*YsumP);
			//Atauimu*rho[i]*min(YFueli, YO2si);
			//Ataui*rhoi*min(YFueli, YO2si);
		if (x[i] < -1.56e-2)  // no combustion in the inlet pipe
		{
			this->wFuel_[i] = 0.0;
		}
	    }
	    const fvPatchList& patches = this->mesh().boundary();
    	    const labelList& own = this->mesh().owner();
	    #include "inertWalls.H" // reaction rate to zero at first cell layer w.r.t walls, scale down second layer
	    //this->wFuel_ = wFuels;
	    /*Info<< minAtauimu << " < Atauimu < " << maxAtauimu << endl;
	    Info<< minAtauik << " < Atauik < " << maxAtauik << endl;
	    Info<< min(this->wFuel_).value() << " < this->wFuel_ < " << max(this->wFuel_).value() << endl;*/
        }
    }
}


template<class CombThermoType, class ThermoType>
bool EDM<CombThermoType, ThermoType>::read()
{
    if (singleStepCombustion<CombThermoType, ThermoType>::read())
    {
        this->coeffs().lookup("C") >> C_ ;
        return true;
    }
    else
    {
        return false;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace combustionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
