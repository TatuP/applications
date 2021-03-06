/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012 OpenFOAM Foundation
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

#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "slicedVolFields.H"
#include "slicedSurfaceFields.H"
#include "SubField.H"
#include "demandDrivenData.H"
#include "fvMeshLduAddressing.H"
#include "emptyPolyPatch.H"
#include "mapPolyMesh.H"
#include "MapFvFields.H"
#include "fvMeshMapper.H"
#include "mapClouds.H"

#include "volPointInterpolation.H"
#include "extendedLeastSquaresVectors.H"
#include "extendedLeastSquaresVectors.H"
#include "leastSquaresVectors.H"
#include "CentredFitData.H"
#include "linearFitPolynomial.H"
#include "quadraticFitPolynomial.H"
#include "quadraticLinearFitPolynomial.H"
//#include "quadraticFitSnGradData.H"
#include "skewCorrectionVectors.H"


#include "centredCECCellToFaceStencilObject.H"
#include "centredCFCCellToFaceStencilObject.H"
#include "centredCPCCellToFaceStencilObject.H"
#include "centredFECCellToFaceStencilObject.H"
#include "upwindCECCellToFaceStencilObject.H"
#include "upwindCFCCellToFaceStencilObject.H"
#include "upwindCPCCellToFaceStencilObject.H"
#include "upwindFECCellToFaceStencilObject.H"

#include "centredCFCFaceToCellStencilObject.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::fvMesh, 0);

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMesh::clearGeomNotOldVol()
{
    slicedVolScalarField::DimensionedInternalField* VPtr =
        static_cast<slicedVolScalarField::DimensionedInternalField*>(VPtr_);
    deleteDemandDrivenData(VPtr);
    VPtr_ = NULL;

    deleteDemandDrivenData(SfPtr_);
    deleteDemandDrivenData(magSfPtr_);
    deleteDemandDrivenData(CPtr_);
    deleteDemandDrivenData(CfPtr_);
}


void Foam::fvMesh::updateGeomNotOldVol()
{
    bool haveV = (VPtr_ != NULL);
    bool haveSf = (SfPtr_ != NULL);
    bool haveMagSf = (magSfPtr_ != NULL);
    bool haveCP = (CPtr_ != NULL);
    bool haveCf = (CfPtr_ != NULL);

    clearGeomNotOldVol();

    // Now recreate the fields
    if (haveV)
    {
        (void)V();
    }
    if (haveSf)
    {
        (void)Sf();
    }
    if (haveMagSf)
    {
        (void)magSf();
    }
    if (haveCP)
    {
        (void)C();
    }
    if (haveCf)
    {
        (void)Cf();
    }
}


void Foam::fvMesh::clearGeom()
{
    clearGeomNotOldVol();

    deleteDemandDrivenData(V0Ptr_);
    deleteDemandDrivenData(V00Ptr_);

    // Mesh motion flux cannot be deleted here because the old-time flux
    // needs to be saved.

    // Things geometry dependent that are not updated.
    volPointInterpolation::Delete(*this);
    extendedLeastSquaresVectors::Delete(*this);
    leastSquaresVectors::Delete(*this);
    CentredFitData<linearFitPolynomial>::Delete(*this);
    CentredFitData<quadraticFitPolynomial>::Delete(*this);
    CentredFitData<quadraticLinearFitPolynomial>::Delete(*this);
    skewCorrectionVectors::Delete(*this);
    //quadraticFitSnGradData::Delete(*this);
}


void Foam::fvMesh::clearAddressing()
{
    deleteDemandDrivenData(lduPtr_);

    // Hack until proper callbacks. Below are all the fvMesh-MeshObjects.

    volPointInterpolation::Delete(*this);
    extendedLeastSquaresVectors::Delete(*this);
    leastSquaresVectors::Delete(*this);
    CentredFitData<linearFitPolynomial>::Delete(*this);
    CentredFitData<quadraticFitPolynomial>::Delete(*this);
    CentredFitData<quadraticLinearFitPolynomial>::Delete(*this);
    skewCorrectionVectors::Delete(*this);
    //quadraticFitSnGradData::Delete(*this);

    centredCECCellToFaceStencilObject::Delete(*this);
    centredCFCCellToFaceStencilObject::Delete(*this);
    centredCPCCellToFaceStencilObject::Delete(*this);
    centredFECCellToFaceStencilObject::Delete(*this);
    // Is this geometry related - cells distorting to upwind direction?
    upwindCECCellToFaceStencilObject::Delete(*this);
    upwindCFCCellToFaceStencilObject::Delete(*this);
    upwindCPCCellToFaceStencilObject::Delete(*this);
    upwindFECCellToFaceStencilObject::Delete(*this);

    centredCFCFaceToCellStencilObject::Delete(*this);
}


void Foam::fvMesh::clearOut()
{
    clearGeom();
    surfaceInterpolation::clearOut();

    clearAddressing();

    // Clear mesh motion flux
    deleteDemandDrivenData(phiPtr_);

    polyMesh::clearOut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMesh::fvMesh(const IOobject& io)
:
    polyMesh(io),
    surfaceInterpolation(*this),
    fvSchemes(static_cast<const objectRegistry&>(*this)),
    fvSolution(static_cast<const objectRegistry&>(*this)),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this, boundaryMesh()),
    lduPtr_(NULL),
    curTimeIndex_(time().timeIndex()),
    VPtr_(NULL),
    V0Ptr_(NULL),
    V00Ptr_(NULL),
    SfPtr_(NULL),
    magSfPtr_(NULL),
    CPtr_(NULL),
    CfPtr_(NULL),
    phiPtr_(NULL)
{
    if (debug)
    {
        Info<< "Constructing fvMesh from IOobject"
            << endl;
    }

    // Check the existance of the cell volumes and read if present
    // and set the storage of V00
    if (isFile(time().timePath()/"V0"))
    {
        V0Ptr_ = new DimensionedField<scalar, volMesh>
        (
            IOobject
            (
                "V0",
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            *this
        );

        V00();
    }

    // Check the existance of the mesh fluxes, read if present and set the
    // mesh to be moving
    if (isFile(time().timePath()/"meshPhi"))
    {
        phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                time().timeName(),
                *this,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            *this
        );

        // The mesh is now considered moving so the old-time cell volumes
        // will be required for the time derivatives so if they haven't been
        // read initialise to the current cell volumes
        if (!V0Ptr_)
        {
            V0Ptr_ = new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "V0",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                V()
            );
        }

        moving(true);
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const Xfer<faceList>& faces,
    const Xfer<labelList>& allOwner,
    const Xfer<labelList>& allNeighbour,
    const bool syncPar
)
:
    polyMesh(io, points, faces, allOwner, allNeighbour, syncPar),
    surfaceInterpolation(*this),
    fvSchemes(static_cast<const objectRegistry&>(*this)),
    fvSolution(static_cast<const objectRegistry&>(*this)),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this),
    lduPtr_(NULL),
    curTimeIndex_(time().timeIndex()),
    VPtr_(NULL),
    V0Ptr_(NULL),
    V00Ptr_(NULL),
    SfPtr_(NULL),
    magSfPtr_(NULL),
    CPtr_(NULL),
    CfPtr_(NULL),
    phiPtr_(NULL)
{
    if (debug)
    {
        Info<< "Constructing fvMesh from components" << endl;
    }
}


Foam::fvMesh::fvMesh
(
    const IOobject& io,
    const Xfer<pointField>& points,
    const Xfer<faceList>& faces,
    const Xfer<cellList>& cells,
    const bool syncPar
)
:
    polyMesh(io, points, faces, cells, syncPar),
    surfaceInterpolation(*this),
    fvSchemes(static_cast<const objectRegistry&>(*this)),
    fvSolution(static_cast<const objectRegistry&>(*this)),
    data(static_cast<const objectRegistry&>(*this)),
    boundary_(*this),
    lduPtr_(NULL),
    curTimeIndex_(time().timeIndex()),
    VPtr_(NULL),
    V0Ptr_(NULL),
    V00Ptr_(NULL),
    SfPtr_(NULL),
    magSfPtr_(NULL),
    CPtr_(NULL),
    CfPtr_(NULL),
    phiPtr_(NULL)
{
    if (debug)
    {
        Info<< "Constructing fvMesh from components" << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fvMesh::~fvMesh()
{
    clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fvMesh::addFvPatches
(
    const List<polyPatch*> & p,
    const bool validBoundary
)
{
    if (boundary().size())
    {
        FatalErrorIn
        (
            "fvMesh::addFvPatches(const List<polyPatch*>&, const bool)"
        )   << " boundary already exists"
            << abort(FatalError);
    }

    // first add polyPatches
    addPatches(p, validBoundary);
    boundary_.addPatches(boundaryMesh());
}


void Foam::fvMesh::removeFvBoundary()
{
    if (debug)
    {
        Info<< "void fvMesh::removeFvBoundary(): "
            << "Removing boundary patches."
            << endl;
    }

    // Remove fvBoundaryMesh data first.
    boundary_.clear();
    boundary_.setSize(0);
    polyMesh::removeBoundary();

    clearOut();
}


Foam::polyMesh::readUpdateState Foam::fvMesh::readUpdate()
{
    if (debug)
    {
        Info<< "polyMesh::readUpdateState fvMesh::readUpdate() : "
            << "Updating fvMesh.  ";
    }

    polyMesh::readUpdateState state = polyMesh::readUpdate();

    if (state == polyMesh::TOPO_PATCH_CHANGE)
    {
        if (debug)
        {
            Info<< "Boundary and topological update" << endl;
        }

        boundary_.readUpdate(boundaryMesh());

        clearOut();

    }
    else if (state == polyMesh::TOPO_CHANGE)
    {
        if (debug)
        {
            Info<< "Topological update" << endl;
        }

        clearOut();
    }
    else if (state == polyMesh::POINTS_MOVED)
    {
        if (debug)
        {
            Info<< "Point motion update" << endl;
        }

        clearGeom();
    }
    else
    {
        if (debug)
        {
            Info<< "No update" << endl;
        }
    }

    return state;
}


const Foam::fvBoundaryMesh& Foam::fvMesh::boundary() const
{
    return boundary_;
}


const Foam::lduAddressing& Foam::fvMesh::lduAddr() const
{
    if (!lduPtr_)
    {
        lduPtr_ = new fvMeshLduAddressing(*this);
    }

    return *lduPtr_;
}


void Foam::fvMesh::mapFields(const mapPolyMesh& meshMap)
{
    // Create a mapper
    const fvMeshMapper mapper(*this, meshMap);

    // Map all the volFields in the objectRegistry
    MapGeometricFields<scalar, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<vector, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<sphericalTensor, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<symmTensor, fvPatchField, fvMeshMapper, volMesh>
    (mapper);
    MapGeometricFields<tensor, fvPatchField, fvMeshMapper, volMesh>
    (mapper);

    // Map all the surfaceFields in the objectRegistry
    MapGeometricFields<scalar, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields<vector, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields<symmTensor, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields<symmTensor, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);
    MapGeometricFields<tensor, fvsPatchField, fvMeshMapper, surfaceMesh>
    (mapper);

    // Map all the clouds in the objectRegistry
    mapClouds(*this, meshMap);


    const labelList& cellMap = meshMap.cellMap();

    // Map the old volume. Just map to new cell labels.
    if (V0Ptr_)
    {
        scalarField& V0 = *V0Ptr_;

        scalarField savedV0(V0);
        V0.setSize(nCells());

        forAll(V0, i)
        {
            if (cellMap[i] > -1)
            {
                V0[i] = savedV0[cellMap[i]];
            }
            else
            {
                V0[i] = 0.0;
            }
        }
    }

    // Map the old-old volume. Just map to new cell labels.
    if (V00Ptr_)
    {
        scalarField& V00 = *V00Ptr_;

        scalarField savedV00(V00);
        V00.setSize(nCells());

        forAll(V00, i)
        {
            if (cellMap[i] > -1)
            {
                V00[i] = savedV00[cellMap[i]];
            }
            else
            {
                V00[i] = 0.0;
            }
        }
    }
}


// Temporary helper function to call move points on
// MeshObjects
template<class Type>
void MeshObjectMovePoints(const Foam::fvMesh& mesh)
{
    if (mesh.thisDb().foundObject<Type>(Type::typeName))
    {
        const_cast<Type&>
        (
            mesh.thisDb().lookupObject<Type>
            (
                Type::typeName
            )
        ).movePoints();
    }
}


Foam::tmp<Foam::scalarField> Foam::fvMesh::movePoints(const pointField& p)
{
    // Grab old time volumes if the time has been incremented
    if (curTimeIndex_ < time().timeIndex())
    {
        if (V00Ptr_ && V0Ptr_)
        {
            *V00Ptr_ = *V0Ptr_;
        }

        if (V0Ptr_)
        {
            *V0Ptr_ = V();
        }
        else
        {
            V0Ptr_ = new DimensionedField<scalar, volMesh>
            (
                IOobject
                (
                    "V0",
                    time().timeName(),
                    *this,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                V()
            );
        }

        curTimeIndex_ = time().timeIndex();
    }


    if (!phiPtr_)
    {
        // Create mesh motion flux
        phiPtr_ = new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                this->time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            *this,
            dimVolume/dimTime
        );
    }
    else
    {
        // Grab old time mesh motion fluxes if the time has been incremented
        if (phiPtr_->timeIndex() != time().timeIndex())
        {
            phiPtr_->oldTime();
        }
    }

    surfaceScalarField& phi = *phiPtr_;

    // Move the polyMesh and set the mesh motion fluxes to the swept-volumes

    scalar rDeltaT = 1.0/time().deltaTValue();

    tmp<scalarField> tsweptVols = polyMesh::movePoints(p);
    scalarField& sweptVols = tsweptVols();

    phi.internalField() = scalarField::subField(sweptVols, nInternalFaces());
    phi.internalField() *= rDeltaT;

    const fvPatchList& patches = boundary();

    forAll(patches, patchI)
    {
        phi.boundaryField()[patchI] = patches[patchI].patchSlice(sweptVols);
        phi.boundaryField()[patchI] *= rDeltaT;
    }

    // Update or delete the local geometric properties as early as possible so
    // they can be used if necessary. These get recreated here instead of
    // demand driven since they might do parallel transfers which can conflict
    // with when they're actually being used.
    // Note that between above "polyMesh::movePoints(p)" and here nothing
    // should use the local geometric properties.
    updateGeomNotOldVol();


    // Update other local data
    boundary_.movePoints();
    surfaceInterpolation::movePoints();


    // Hack until proper callbacks. Below are all the fvMesh MeshObjects with a
    // movePoints function.
    MeshObjectMovePoints<volPointInterpolation>(*this);
    MeshObjectMovePoints<extendedLeastSquaresVectors>(*this);
    MeshObjectMovePoints<leastSquaresVectors>(*this);
    MeshObjectMovePoints<CentredFitData<linearFitPolynomial> >(*this);
    MeshObjectMovePoints<CentredFitData<quadraticFitPolynomial> >(*this);
    MeshObjectMovePoints<CentredFitData<quadraticLinearFitPolynomial> >(*this);
    MeshObjectMovePoints<skewCorrectionVectors>(*this);
    //MeshObjectMovePoints<quadraticFitSnGradData>(*this);

    return tsweptVols;
}


void Foam::fvMesh::updateMesh(const mapPolyMesh& mpm)
{
    // Update polyMesh. This needs to keep volume existent!
    polyMesh::updateMesh(mpm);

    // Clear the sliced fields
    clearGeomNotOldVol();

    // Map all fields using current (i.e. not yet mapped) volume
    mapFields(mpm);

    // Clear the current volume and other geometry factors
    surfaceInterpolation::clearOut();

    clearAddressing();

    // handleMorph() should also clear out the surfaceInterpolation.
    // This is a temporary solution
    surfaceInterpolation::movePoints();
}


bool Foam::fvMesh::writeObjects
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    return polyMesh::writeObject(fmt, ver, cmp);
}


//- Write mesh using IO settings from the time
bool Foam::fvMesh::write() const
{
    return polyMesh::write();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

bool Foam::fvMesh::operator!=(const fvMesh& bm) const
{
    return &bm != this;
}


bool Foam::fvMesh::operator==(const fvMesh& bm) const
{
    return &bm == this;
}


// ************************************************************************* //
