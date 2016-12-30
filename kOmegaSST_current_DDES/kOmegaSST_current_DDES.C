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

#include "kOmegaSST_current_DDES.H"
#include "addToRunTimeSelectionTable.H"

#include "backwardsCompatibilityWallFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kOmegaSST_current_DDES, 0);
addToRunTimeSelectionTable(RASModel, kOmegaSST_current_DDES, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> kOmegaSST_current_DDES::F1(const volScalarField& CDkOmega) const
{
    tmp<volScalarField> CDkOmegaPlus = max
    (
        CDkOmega,
        dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
                scalar(500)*nu()/(sqr(y_)*omega_)
            ),
            (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

tmp<volScalarField> kOmegaSST_current_DDES::F2() const
{
    tmp<volScalarField> arg2 = min
    (
        max
        (
            (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
            scalar(500)*nu()/(sqr(y_)*omega_)
        ),
        scalar(100)
    );

    return tanh(sqr(arg2));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaSST_current_DDES::kOmegaSST_current_DDES
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85034
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.85616
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            0.5532
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.4403
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    CDES1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDES1",
            coeffDict_,
            0.78
        )
    ),
    CDES2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDES2",
            coeffDict_,
            0.12
        )
    ),
    Cd1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd1",
            coeffDict_,
            8.0
        )
    ),
    Cd2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cd2",
            coeffDict_,
            3.0
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),

    y_(mesh_),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateK("k", mesh_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateOmega("omega", mesh_)
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        autoCreateNut("nut", mesh_)
    ),
    nut_RANS_
    (
        IOobject
        (
            "nut_RANS",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
        dimensionedScalar("nut_RANS", dimLength*dimVelocity, 0.0)
    ),
    lDDES_
    (
        IOobject
        (
            "lDDES",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        dimensionedScalar("lDDES", dimLength, 0.0)
    ),
    lRANS_
    (
        IOobject
        (
            "lRANS",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        dimensionedScalar("lRANS", dimLength, 0.0)
    ),
   lLESRANS_
    (
        IOobject
        (
            "lLESRANS",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        dimensionedScalar("lLESRANS", dimLength/dimLength, 0.0)
    ),
    lLES_
    (
        IOobject
        (
            "lLES",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        dimensionedScalar("lLES", dimLength, 0.0)
    ),
    lzero_
    (
        IOobject
        (
            "lzero",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
        dimensionedScalar("lzero", dimLength, 0.0)
    ),
    f2s_
    (
        IOobject
        (
            "f2s",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
        dimensionedScalar("f2s", dimVelocity/dimLength, 0.0)
    ),
    fd_
    (
        IOobject
        (
            "fd",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh_,
        dimensionedScalar("fd", dimVelocity/dimVelocity, 0.0)
    ),
    rd_
    (
        IOobject
        (
            "rd",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
        dimensionedScalar("rd", dimVelocity/dimVelocity, 0.0)
    ),
    hmax_
    (
        IOobject
        (
            "hmax",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
        dimensionedScalar("hmax", dimLength, 0.0)
    ),
    delta_
    (
        IOobject
        (
            "delta",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
	mesh_,
        dimensionedScalar("delta", dimLength, 0.0)
    ),
    absS_
    (
        IOobject
        (
            "S",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
    mesh_,
        dimensionedScalar("S", dimVelocity/dimLength, 0.0)
    )
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);

    update_nut();
/*    nut_ =
    (
        a1_*k_
      / max
        (
            a1_*omega_,
            F2()*sqrt(2.0)*mag(symm(fvc::grad(U_)))
        )
    );
    nut_.correctBoundaryConditions();*/

    calchmax();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kOmegaSST_current_DDES::R() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "R",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            ((2.0/3.0)*I)*k_ - nut_*twoSymm(fvc::grad(U_)),
            k_.boundaryField().types()
        )
    );
}


tmp<volSymmTensorField> kOmegaSST_current_DDES::devReff() const
{
    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            IOobject
            (
                "devRhoReff",
                runTime_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
           -nuEff()*dev(twoSymm(fvc::grad(U_)))
        )
    );
}


tmp<fvVectorMatrix> kOmegaSST_current_DDES::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}

tmp<fvVectorMatrix> kOmegaSST_current_DDES::divDevRhoReff
(
    const volScalarField& rho,
    volVectorField& U
) const
{
    volScalarField muEff("muEff", rho*nuEff());

    return
    (
      - fvm::laplacian(muEff, U)
      - fvc::div(muEff*dev(T(fvc::grad(U))))
    );
}

bool kOmegaSST_current_DDES::read()
{
    if (RASModel::read())
    {
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kOmegaSST_current_DDES::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    if (mesh_.changing())
    {
        y_.correct();
    }

    const volScalarField S2(2*magSqr(symm(fvc::grad(U_))));
    const volScalarField O2(2*magSqr(skew(fvc::grad(U_))));
//    strain_ = sqrt(S2);
//    vorticity_ = sqrt(O2);
//    q_crit_ = scalar(-0.5)*(S2 - O2);

    rd_ = (nut_RANS_ + nu())/(sqr(kappa_) * sqr(y_) * sqrt(scalar(0.5)*(S2 + O2)));
    fd_ = scalar(1) - tanh(pow(Cd1_*rd_, Cd2_));
    update_nut();

    volScalarField G("kOmegaSST_current_DDES:G", nut_*S2);

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    const volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1(this->F1(CDkOmega));

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_), omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)*S2
      - fvm::Sp(beta(F1)*omega_, omega_)
      - fvm::SuSp
        (
            (F1 - scalar(1))*CDkOmega/omega_,
            omega_
        )
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G, c1_*betaStar_*k_*omega_)
      - fvm::Sp(betaStar_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);


    // Re-calculate viscosity
//    nut_ = a1_*k_/max(a1_*omega_, F2()*sqrt(S2));
    Info << "min k = " << min(k_)<< "  max k = " << max(k_) << endl;
    Info << "min omega = " << min(omega_)<< "  max omega = " << max(omega_) << endl;
    Info << "min nut = " << min(nut_)<< "  max nut = " << max(nut_) << endl;

}

void kOmegaSST_current_DDES::update_nut()
{
    const volScalarField CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    const volScalarField F1(this->F1(CDkOmega));

    lRANS_ = sqrt(a1_*k_)/max(sqrt(a1_*omega_*omega_), sqrt(omega_*F2()*sqrt(2.0)*mag(symm(fvc::grad(U_)))));
    f2s_ = F2()*sqrt(2.0)*mag(symm(fvc::grad(U_)));
//    lRANS_ = sqrt(k_)/omega_;
    lLES_ = CDES2_*(fd_*delta_ + (scalar(1) - fd_)*hmax_);
    lDDES_ = lRANS_ - fd_*max(lzero_, lRANS_ - lLES_);
    lLESRANS_ = max(lzero_, lRANS_ - lLES_)/(lRANS_ - lLES_); 
    absS_ = sqrt(2*magSqr(symm(fvc::grad(U_))));
    nut_ = lDDES_ * lDDES_ * omega_;
    nut_RANS_ =
    (
        a1_*k_
      / max
        (
            a1_*omega_,
            F2()*sqrt(2.0)*mag(symm(fvc::grad(U_)))
        )
    );
    nut_.correctBoundaryConditions();
}

void kOmegaSST_current_DDES::calchmax()
{
        label nD = mesh_.nGeometricD();
        if(nD==2)
        {
                WarningIn("kOmegaDES::calcDelta()")
                        << "Case is 2D, LES is not strictly applicable\n"
                        << endl;
        }

        const cellList& cells = mesh_.cells();
        forAll(cells,cellI)
        {
                scalar deltaMaxTmp = 0.0;
                const labelList& cFaces = mesh_.cells()[cellI];
                const point& centrevector = mesh_.cellCentres()[cellI];
                forAll(cFaces, cFaceI)
                {
                        label faceI = cFaces[cFaceI];
                        const point& facevector = mesh_.faceCentres()[faceI];
                        scalar tmp = 2*mag(facevector - centrevector);
                        deltaMaxTmp = max(deltaMaxTmp, tmp);
                }
                hmax_[cellI] = deltaMaxTmp;
        }
        delta_.internalField() = pow(mesh_.V(), 1.0/3.0);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
