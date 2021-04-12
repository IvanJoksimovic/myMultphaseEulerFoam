/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::diameterModels::IATE

    Description
        IATE (Interfacial Area Transport Equation) bubble diameter model.

        Solves for the interfacial curvature per unit volume of the phase rather
        than interfacial area per unit volume to avoid stability issues relating to
        the consistency requirements between the phase fraction and interfacial area
        per unit volume.  In every other respect this model is as presented in the
        paper:

        Reference:
        \verbatim
            Ishii, M., Kim, S., & Kelly, J. (2005).
            Development of interfacial area transport equation.
            Nuclear Engineering and Technology, 37(6), 525-536.
        \endverbatim

    SourceFiles
        IATE.C

\*---------------------------------------------------------------------------*/

#ifndef diameterModels_IATE_H
#define diameterModels_IATE_H

#include "diameterModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace diameterModels
{

// Forward declaration of classes
class IATEsource;

/*---------------------------------------------------------------------------*\
                           Class IATE Declaration
\*---------------------------------------------------------------------------*/

class IATE
:
    public diameterModel
{
    // Private Data

    //- Interfacial curvature (alpha*interfacial area)
    volScalarField kappai_;

    //- Maximum diameter used for stabilisation in the limit kappai->0
    dimensionedScalar dMax_;

    //- Minimum diameter used for stabilisation in the limit kappai->inf
    dimensionedScalar dMin_;

    //- Residual phase fraction
    dimensionedScalar residualAlpha_;

    //- The Sauter-mean diameter of the phase
    //volScalarField& d_;

    //- IATE sources
    //PtrList<IATEsource> sources_;




public:

    //- Runtime type information
    TypeName("IATE");


    // Constructors

        //- Construct from components
        IATE
        (
            const dictionary& dict,
            const phaseModel& phase
        );


    //- Destructor
    virtual ~IATE() = default;


    // Member Functions

    tmp<volScalarField> d() const;

    //- Correct the diameter field
    void correct() const;



};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace diameterModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
