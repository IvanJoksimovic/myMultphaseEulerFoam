/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "IATE.H"
#include "IATEsource.H"
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmSup.H"
#include "fvcDdt.H"
#include "fvcDiv.H"
#include "fvcAverage.H"
#include "fvOptions.H"
#include "mathematicalConstants.H"
#include "fundamentalConstants.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{
    defineTypeNameAndDebug(IATE, 0);

    addToRunTimeSelectionTable
    (
        diameterModel,
        IATE,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diameterModels::IATE::IATE
(
    const dictionary& dict,
    const phaseModel& phase
)
:
    diameterModel(dict, phase),
    kappai_
    (
        IOobject
        (
            IOobject::groupName("kappai", phase.name()),
            phase.time().timeName(),
            phase.mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        phase.mesh()
    ),
    dMax_("dMax", dimLength, dict),
    dMin_("dMax", dimLength, dict),
    residualAlpha_("residualAlpha", dimless, dict)//,
    //    sources_(diameterProperties.lookup("sources"), IATEsource::iNew(*this))

{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
Foam::tmp<Foam::volScalarField> Foam::diameterModels::IATE::d() const
{
    correct();

    volScalarField d_(max(6/max(kappai_, 6/dMax_), dMin_)) ;

    Info << "Global dsmr = " << gAverage(d_) << endl;

    Info << phase().db().sortedNames(		) << endl;

    return max(6/max(kappai_, 6/dMax_), dMin_);
}

void Foam::diameterModels::IATE::correct() const
{
    volScalarField alphaAv
    (
        max
        (
            0.5*fvc::average(phase() + phase().oldTime()),
            residualAlpha_
        )
    );

    // Initialise the accumulated source term to the dilatation effect
    // Unnecesary for incompressible flows, but we need to initialise R somewhere
    fvScalarMatrix R
    (
       -fvm::SuSp
        (

            ((1.0/3.0)/alphaAv)
           *(
                (
                    fvc::ddt(phase())
                  + fvc::div(phase().alphaPhi())
                )
              - (
                  fvc::ddt(phase())
                + fvc::div(phase().alphaPhi())
                )
            ),
            kappai_
        )
    );

    // Accumulate the run-time selectable sources
    /*
    forAll(sources_, j)
    {
        R += sources_[j].R(alphaAv, kappai_);
    }
    */

    fv::options& fvOptions(fv::options::New(phase().mesh()));

    // Construct the interfacial curvature equation
    fvScalarMatrix kappaiEqn
    (
        fvm::ddt(kappai_) + fvm::div(phase().phi(), kappai_)
      - fvm::Sp(fvc::div(phase().phi()), kappai_)
     ==
      R
      //+ fvOptions(kappai_)
    );

    kappaiEqn.relax();

    fvOptions.constrain(kappaiEqn);

    kappaiEqn.solve();

    // Update the Sauter-mean diameter
    //d_ = dsm();
}




// ************************************************************************* //
