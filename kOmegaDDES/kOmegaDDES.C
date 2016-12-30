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

#include "kOmegaDDES.H"
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

defineTypeNameAndDebug(kOmegaDDES, 0);
addToRunTimeSelectionTable(RASModel, kOmegaDDES, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kOmegaDDES::kOmegaDDES
(
    const volVectorField& U,
    const surfaceScalarField& phi,
    transportModel& transport,
    const word& turbulenceModelName,
    const word& modelName
)
:
    RASModel(modelName, U, phi, transport, turbulenceModelName),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
//            "betaStar",
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    beta_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
//            "beta",
            "C_omega2",
            coeffDict_,
            0.072
        )
    ),
    alpha_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
//            "alpha",
            "C_omega1",
            coeffDict_,
            0.52
        )
    ),
    alphaK_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega",
            coeffDict_,
            0.5
        )
    ),

    Cu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cu",
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
    CDES_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "CDES",
            coeffDict_,
            0.12
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
    lDDES_
    (
        IOobject
        (
            "lDDES",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
        ),
	mesh_,
	dimensionedScalar("lRANS", dimLength, 0.0)
    ),
    lLES_
    (
        IOobject
        (
            "lLES",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
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

    y_(mesh_)
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);

    calcdelta();
    update_nut();

    printCoeffs();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volSymmTensorField> kOmegaDDES::R() const
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


tmp<volSymmTensorField> kOmegaDDES::devReff() const
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


tmp<fvVectorMatrix> kOmegaDDES::divDevReff(volVectorField& U) const
{
    return
    (
      - fvm::laplacian(nuEff(), U)
      - fvc::div(nuEff()*dev(T(fvc::grad(U))))
    );
}

tmp<fvVectorMatrix> kOmegaDDES::divDevRhoReff
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

bool kOmegaDDES::read()
{
    if (RASModel::read())
    {
        Cmu_.readIfPresent(coeffDict());
        beta_.readIfPresent(coeffDict());
        alphaK_.readIfPresent(coeffDict());
        alphaOmega_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void kOmegaDDES::update_nut()
{
    lRANS_ = sqrt(k_)/(omega_);
    lLES_ = CDES_ * (fd_*delta_ + (scalar(1) - fd_)*hmax_);
    lDDES_ = lRANS_ - fd_ * max(lzero_, lRANS_ - lLES_);
    nut_ = lDDES_ * lDDES_ * omega_;
    nut_.correctBoundaryConditions();
}

void kOmegaDDES::correct()
{
    RASModel::correct();

    if (!turbulence_)
    {
        return;
    }

    const volScalarField S2(2*magSqr(symm(fvc::grad(U_))));
    const volScalarField O2(2*magSqr(skew(fvc::grad(U_))));
//    rd_ = (nut_ + nu())/(sqr(kappa_) * sqr(y_) * sqrt(scalar(0.5)*(S2 + O2)));
//--------------k/omega instead of nut for rd---------------------------
    rd_ = ((k_/omega_) + nu())/(sqr(kappa_) * sqr(y_) * sqrt(scalar(0.5)*(S2 + O2)));
    fd_ = scalar(1) - tanh(pow(Cd1_*rd_, Cd2_));
    // Re-calculate viscosity
    update_nut();

    volScalarField G("kOmegaDDES:G", nut_*2*magSqr(symm(fvc::grad(U_))));

    // Update omega and G at the wall
    omega_.boundaryField().updateCoeffs();

    // Turbulence specific dissipation rate equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::Sp(fvc::div(phi_), omega_)
      - fvm::laplacian(DomegaEff(), omega_)
     ==
//        alpha_*G*omega_/k_
        alpha_*2*magSqr(symm(fvc::grad(U_)))
      - fvm::Sp(beta_*omega_, omega_)
    );

    omegaEqn().relax();

    omegaEqn().boundaryManipulate(omega_.boundaryField());

    solve(omegaEqn);
    bound(omega_, omegaMin_);
    omega_ = max(omega_, scalar(0.02)*sqrt(2*magSqr(symm(fvc::grad(U_)))));

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G
      - fvm::Sp(Cmu_*omega_, k_)
    );

    kEqn().relax();
    solve(kEqn);
    bound(k_, kMin_);

  //  Info << "min k = " << min(k_)<< "  max k = " << max(k_) << endl;
  //  Info << "min omega = " << min(omega_)<< "  max omega = " << max(omega_) << endl;
  //  Info << "min hmax = " << min(hmax_)<< "  max hmax = " << max(hmax_) << endl;
  //  Info << "min nut = " << min(nut_)<< "  max nut = " << max(nut_) << endl;

}

void kOmegaDDES::calcdelta()
{
	label nD = mesh_.nGeometricD();
	if(nD==2)
	{
		WarningIn("kOmegaDDES::calcDelta()")
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
