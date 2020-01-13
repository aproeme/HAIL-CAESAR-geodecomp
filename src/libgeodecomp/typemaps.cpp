#include "typemaps.h"
#include <algorithm>
#include <stdexcept>

MPI_Datatype MPI_CELL;


// Member Specification, holds all relevant information for a given member.
class MemberSpec
{
public:
    MemberSpec(MPI_Aint address, MPI_Datatype type, int length) :
        address(address),
        type(type),
        length(length)
    {}

    MPI_Aint address;
    MPI_Datatype type;
    int length;
};

bool addressLower(MemberSpec a, MemberSpec b)
{
    return a.address < b.address;
}

MPI_Datatype
Typemaps::generateMapCell() {
    char fakeObject[sizeof(Cell)];
    Cell *obj = (Cell*)fakeObject;

    const int count = 5;
    int lengths[count];

    // sort addresses in ascending order
    MemberSpec rawSpecs[] = {
        MemberSpec(getAddress(&obj->celltype), lookup<int >(), 1),
        MemberSpec(getAddress(&obj->elevation), lookup<double >(), 1),
        MemberSpec(getAddress(&obj->qx), lookup<double >(), 1),
        MemberSpec(getAddress(&obj->qy), lookup<double >(), 1),
        MemberSpec(getAddress(&obj->water_depth), lookup<double >(), 1)
    };
    std::sort(rawSpecs, rawSpecs + count, addressLower);

    // split addresses from member types
    MPI_Aint displacements[count];
    MPI_Datatype memberTypes[count];
    for (int i = 0; i < count; i++) {
        displacements[i] = rawSpecs[i].address;
        memberTypes[i] = rawSpecs[i].type;
        lengths[i] = rawSpecs[i].length;
    }

    // transform absolute addresses into offsets
    for (int i = count-1; i > 0; i--) {
        displacements[i] -= displacements[0];
    }
    displacements[0] = 0;

    // create datatype
    MPI_Datatype objType;
    MPI_Type_create_struct(count, lengths, displacements, memberTypes, &objType);
    MPI_Type_commit(&objType);

    return objType;
}


void Typemaps::initializeMaps()
{
    if (mapsCreated) {
        throw std::logic_error("Typemaps already initialized, duplicate initialization would leak memory");
    }

    if (sizeof(std::size_t) != sizeof(unsigned long)) {
        throw std::logic_error("MPI_UNSIGNED_LONG not suited for communication of std::size_t, needs to be redefined");
    }

    int mpiInitState = 0;
    MPI_Initialized(&mpiInitState);
    if (!mpiInitState) {
        throw std::logic_error("MPI needs to be initialized prior to setting up Typemaps");
    }

    MPI_CELL = generateMapCell();

    mapsCreated = true;
}

bool Typemaps::mapsCreated = false;

