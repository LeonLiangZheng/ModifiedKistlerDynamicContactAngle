#include "modifiedKistlerAlphaContactAngleFvPatchScalarField.H"

#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"
#include <cmath>

namespace Foam
{

static const scalar mkSmallThetaRad = 1.0e-12;
static const scalar mkMaxBracket = 1.0e12;
static const label  mkMaxBisectionIter = 200;

scalar modifiedKistlerAlphaContactAngleFvPatchScalarField::degToRad
(
    const scalar deg
)
{
    return deg*constant::mathematical::pi/180.0;
}

scalar modifiedKistlerAlphaContactAngleFvPatchScalarField::radToDeg
(
    const scalar rad
)
{
    return rad*180.0/constant::mathematical::pi;
}

scalar modifiedKistlerAlphaContactAngleFvPatchScalarField::clampDeg
(
    const scalar thetaDeg
)
{
    return max(scalar(0), min(scalar(180), thetaDeg));
}

modifiedKistlerAlphaContactAngleFvPatchScalarField::
modifiedKistlerAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF),
    thetaE_(0.0),
    thetaA_(0.0),
    thetaR_(0.0),
    muEffective_(0.0),
    sigma_(1.0)
{}

modifiedKistlerAlphaContactAngleFvPatchScalarField::
modifiedKistlerAlphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(p, iF, dict),
    thetaE_(dict.get<scalar>("thetaE")),
    thetaA_(dict.get<scalar>("thetaA")),
    thetaR_(dict.get<scalar>("thetaR")),
    muEffective_(dict.get<scalar>("muEffective")),
    sigma_(dict.get<scalar>("sigma"))
{
    evaluate();
}

modifiedKistlerAlphaContactAngleFvPatchScalarField::
modifiedKistlerAlphaContactAngleFvPatchScalarField
(
    const modifiedKistlerAlphaContactAngleFvPatchScalarField& rhs,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(rhs, p, iF, mapper),
    thetaE_(rhs.thetaE_),
    thetaA_(rhs.thetaA_),
    thetaR_(rhs.thetaR_),
    muEffective_(rhs.muEffective_),
    sigma_(rhs.sigma_)
{}

modifiedKistlerAlphaContactAngleFvPatchScalarField::
modifiedKistlerAlphaContactAngleFvPatchScalarField
(
    const modifiedKistlerAlphaContactAngleFvPatchScalarField& rhs
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(rhs),
    thetaE_(rhs.thetaE_),
    thetaA_(rhs.thetaA_),
    thetaR_(rhs.thetaR_),
    muEffective_(rhs.muEffective_),
    sigma_(rhs.sigma_)
{}

modifiedKistlerAlphaContactAngleFvPatchScalarField::
modifiedKistlerAlphaContactAngleFvPatchScalarField
(
    const modifiedKistlerAlphaContactAngleFvPatchScalarField& rhs,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleTwoPhaseFvPatchScalarField(rhs, iF),
    thetaE_(rhs.thetaE_),
    thetaA_(rhs.thetaA_),
    thetaR_(rhs.thetaR_),
    muEffective_(rhs.muEffective_),
    sigma_(rhs.sigma_)
{}

scalar modifiedKistlerAlphaContactAngleFvPatchScalarField::hoffmanFunctionRad
(
    const scalar x
) const
{
    const scalar xx = max(scalar(0), x);
    const scalar inner =
        5.16*pow(xx/(1.0 + 1.31*pow(xx, 0.99)), 0.706);

    const scalar arg = 1.0 - 2.0*std::tanh(inner);
    const scalar boundedArg = max(scalar(-1), min(scalar(1), arg));

    return std::acos(boundedArg);
}

scalar modifiedKistlerAlphaContactAngleFvPatchScalarField::inverseHoffmanRad
(
    const scalar thetaRad
) const
{
    const scalar theta =
        max(mkSmallThetaRad, min(constant::mathematical::pi - mkSmallThetaRad, thetaRad));

    scalar lo = 0.0;
    scalar hi = 1.0;

    while (hoffmanFunctionRad(hi) < theta && hi < mkMaxBracket)
    {
        hi *= 2.0;
    }

    for (label iter = 0; iter < mkMaxBisectionIter; ++iter)
    {
        const scalar mid = 0.5*(lo + hi);
        const scalar fmid = hoffmanFunctionRad(mid);

        if (mag(fmid - theta) < 1.0e-12)
        {
            return mid;
        }

        if (fmid < theta)
        {
            lo = mid;
        }
        else
        {
            hi = mid;
        }
    }

    return 0.5*(lo + hi);
}

scalar modifiedKistlerAlphaContactAngleFvPatchScalarField::modifiedKistlerThetaDeg
(
    const scalar ucl
) const
{
    if (sigma_ <= ROOTVSMALL || muEffective_ < 0)
    {
        return clampDeg(thetaE_);
    }

    const scalar thetaEqRad = degToRad(clampDeg(thetaE_));
    const scalar invEq = inverseHoffmanRad(thetaEqRad);
    const scalar Ca = muEffective_*mag(ucl)/sigma_;

    scalar thetaDynDeg = thetaE_;

    if (ucl >= 0.0)
    {
        const scalar thetaDynRad = hoffmanFunctionRad(Ca + invEq);
        thetaDynDeg = min(thetaA_, radToDeg(thetaDynRad));
    }
    else
    {
        const scalar thetaDynRad = 2.0*thetaEqRad - hoffmanFunctionRad(Ca + invEq);
        thetaDynDeg = max(thetaR_, radToDeg(thetaDynRad));
    }

    return clampDeg(thetaDynDeg);
}

tmp<scalarField>
modifiedKistlerAlphaContactAngleFvPatchScalarField::theta
(
    const fvPatchVectorField& Up,
    const fvsPatchVectorField& nHat
) const
{
    const vectorField nf(patch().nf());

    // Same contact-line velocity extraction pattern as the built-in
    // dynamicAlphaContactAngle: wall-parallel velocity projected onto the
    // wall-parallel interface-normal direction.
    vectorField Uwall(Up.patchInternalField() - Up);
    Uwall -= (nf & Uwall)*nf;

    vectorField nWall(nHat - (nf & nHat)*nf);
    nWall /= (mag(nWall) + ROOTVSMALL);

    const scalarField ucl(nWall & Uwall);

    tmp<scalarField> ttheta(new scalarField(size(), thetaE_));
    scalarField& thetaField = ttheta.ref();

    forAll(thetaField, facei)
    {
        // If you find advancing/receding reversed in your case due to patch
        // orientation, change ucl[facei] to -ucl[facei] here.
        thetaField[facei] = modifiedKistlerThetaDeg(ucl[facei]);
    }

    return ttheta;
}

void modifiedKistlerAlphaContactAngleFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleTwoPhaseFvPatchScalarField::write(os);
    os.writeKeyword("thetaE") << thetaE_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;
    os.writeKeyword("muEffective") << muEffective_ << token::END_STATEMENT << nl;
    os.writeKeyword("sigma") << sigma_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}

makePatchTypeField
(
    fvPatchScalarField,
    modifiedKistlerAlphaContactAngleFvPatchScalarField
);

} // End namespace Foam
