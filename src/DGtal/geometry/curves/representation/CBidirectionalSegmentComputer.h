/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

#pragma once

/**
 * @file CBidirectionalSegmentComputer.h
 * @author Tristan Roussillon (\c tristan.roussillon@liris.cnrs.fr )
 * Laboratoire d'InfoRmatique en Image et Systèmes d'information - LIRIS (CNRS, UMR 5205), CNRS, France
 *
 * @date 2011/08/31
 *
 * Header file for concept CBidirectionalSegmentComputer.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CBidirectionalSegmentComputer_RECURSES)
#error Recursive header files inclusion detected in CBidirectionalSegmentComputer.h
#else // defined(CBidirectionalSegmentComputer_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CBidirectionalSegmentComputer_RECURSES

#if !defined CBidirectionalSegmentComputer_h
/** Prevents repeated inclusion of headers. */
#define CBidirectionalSegmentComputer_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/geometry/curves/representation/CForwardSegmentComputer.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class CBidirectionalSegmentComputer
  /**
Description of \b concept '\b CBidirectionalSegmentComputer' <p>
     @ingroup Concepts
     @brief Aim: Defines the concept describing a bidirectional segment computer,  
    ie. a model of CSegment that can extend itself in the two possible directions. 
     
 ### Refinement of CForwardSegmentComputer 
    
 ### Associated types : the same as CForwardSegmentComputer
  
 ### Notation
     - \t X : A type that is a model of CBidirectionalSegmentComputer
     - \t x : object of type X
  
 ### Definitions
    
 ### Valid expressions and 
     <table> 
      <tr> 
        <td class=CName> \b Name </td> 
        <td class=CExpression> \b Expression </td>
        <td class=CRequirements> \b Type requirements </td> 
        <td class=CReturnType> \b Return type </td>
        <td class=CPrecondition> \b Precondition </td> 
        <td class=CSemantics> \b Semantics </td> 
        <td class=CPostCondition> \b Postcondition </td> 
        <td class=CComplexity> \b Complexity </td>
      </tr>
      <tr> 
        <td class=CName> extension test  </td> 
        <td class=CExpression> x.isExtendableBackward()     </td>
        <td class=CRequirements>    </td> 
        <td class=CReturnType> bool     </td>
        <td class=CPrecondition>    </td> 
        <td class=CSemantics> check wether x can be extended to --x.begin() or not </td> 
        <td class=CPostCondition>       </td> 
        <td class=CComplexity>      </td>
      </tr>
      <tr> 
        <td class=CName> extension </td> 
        <td class=CExpression> x.extendBackward()     </td>
        <td class=CRequirements>    </td> 
        <td class=CReturnType> bool     </td>
        <td class=CPrecondition>    </td> 
        <td class=CSemantics> check wether x can be extended to --x.begin() or not, extend if true </td> 
        <td class=CPostCondition>       </td> 
        <td class=CComplexity>     </td>
      </tr>
     </table>
    
 ### Invariants###
    
 ### Models###

   ArithmeticalDSS3d, GeometricalDSS, GeometricalDCA    

 ### Notes###

@tparam T the type that should be a model of CBidirectionalSegmentComputer.
   */
  template <typename T> 
  struct CBidirectionalSegmentComputer : CForwardSegmentComputer<T>
  {
    // ----------------------- Concept checks ------------------------------
  public:
    // Methods
    BOOST_CONCEPT_USAGE( CBidirectionalSegmentComputer )
    {
      ConceptUtils::sameType( myB, myX.isExtendableBackward() );
      ConceptUtils::sameType( myB, myX.extendBackward() );
    }
    // ------------------------- Private Datas --------------------------------
  private:
    T myX; // only if T is default constructible.
    bool myB; 
  
    // ------------------------- Internals ------------------------------------
  private:
    
  }; // end of concept CBidirectionalSegmentComputer
  
} // namespace DGtal

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CBidirectionalSegmentComputer_h

#undef CBidirectionalSegmentComputer_RECURSES
#endif // else defined(CBidirectionalSegmentComputer_RECURSES)
