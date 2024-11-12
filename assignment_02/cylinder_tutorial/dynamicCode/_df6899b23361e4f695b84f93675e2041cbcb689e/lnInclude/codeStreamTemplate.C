/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) YEAR AUTHOR, AFFILIATION
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "pointField.H"
#include "tensor.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 44 "/work3/s232439/s232439-v2312/run/assignment_02/cylinder/system/blockMeshDict/#codeStream"
#include "pointField.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C" void codeStream_df6899b23361e4f695b84f93675e2041cbcb689e(Foam::Ostream& os, const Foam::dictionary& dict)
{
//{{{ begin code
    #line 49 "/work3/s232439/s232439-v2312/run/assignment_02/cylinder/system/blockMeshDict/#codeStream"
pointField points
        ({
            /* 0*/ { 0.5, 0, -0.5 },
            /* 1*/ { 1, 0, -0.5 },
            /* 2*/ { 2, 0, -0.5 },
            /* 3*/ { 2, 0.7071067811865476, -0.5 },
            /* 4*/ { 0.7071067811865476, 0.7071067811865476, -0.5 },
            /* 5*/ { 0.3535533905932738, 0.3535533905932738, -0.5 },
            /* 6*/ { 2, 2, -0.5 },
            /* 7*/ { 0.7071067811865476, 2, -0.5 },
            /* 8*/ { 0, 2, -0.5 },
            /* 9*/ { 0, 1, -0.5 },
            /*10*/ { 0, 0.5, -0.5 },
            /*11*/ { -0.5, 0, -0.5 },
            /*12*/ { -1, 0, -0.5 },
            /*13*/ { -2, 0, -0.5 },
            /*14*/ { -2, 0.7071067811865476, -0.5 },
            /*15*/ { -0.7071067811865476, 0.7071067811865476, -0.5 },
            /*16*/ { -0.3535533905932738, 0.3535533905932738, -0.5 },
            /*17*/ { -2, 2, -0.5 },
            /*18*/ { -0.7071067811865476, 2, -0.5 }
        });

        // Duplicate z points for zmax
        const label sz = points.size();
        points.resize(2*sz);
        for (label i = 0; i < sz; ++i)
        {
            const point& pt = points[i];
            points[i + sz] = point(pt.x(), pt.y(), 0.5);
        }

        os  << points;
//}}} end code
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

