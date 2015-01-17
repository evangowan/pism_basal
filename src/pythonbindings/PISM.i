// Copyright (C) 2011, 2012, 2013, 2014, 2015 David Maxwell and Constantine Khroulev
//
// This file is part of PISM.
//
// PISM is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// PISM is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License
// along with PISM; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

%module(directors="1") cpp

/* Don't warn about nested classes (325) and methods that can't be
 * wrapped under the present name, i.e. methods that require %rename.
 */
#pragma SWIG nowarn=325,503

%{
// The material in this section is included verbatim in the C++ source code generated by SWIG.
// The necessary header files required to compile must be included.
// This list is NOT the whole set of headers being wrapped; it is just the list of includes that 
// draws in all the other needed includes as well. See the end of this file for the list
// of PISM headers being wrapped.

#include "PISMUnits.hh"
#include "pism_python.hh"

#include "Mask.hh"
#include "basal_resistance.hh"
#include "enthalpyConverter.hh"
#include "PISMMohrCoulombYieldStress.hh"
#include "error_handling.hh"
%}

// Include petsc4py.i so that we get support for automatic handling of PetscErrorCode return values
%include "petsc4py/petsc4py.i"

%include "pism_exception.i"

/* PISM still uses PetscErrorCode in some places. Convert it to integers. */
%typemap(typecheck,precedence=SWIG_TYPECHECK_INTEGER) PetscErrorCode {
   $1 = PyInt_Check($input) ? 1 : 0;
}
%typemap(directorout) PetscErrorCode %{ $result =  PyInt_AsLong($input); %}

// Automatic conversions between std::string and python string arguments and return values
%include std_string.i
// Conversions between python lists and certain std::vector's
%include std_vector.i
%include std_set.i

#ifdef PISM_USE_TR1
#define SWIG_SHARED_PTR_SUBNAMESPACE tr1
#endif
%include <std_shared_ptr.i>

%template(IntVector) std::vector<int>;
%template(DoubleVector) std::vector<double>;
%template(StringVector) std::vector<std::string>;
%template(StringSet) std::set<std::string>;

// Why did I include this?
%include "cstring.i"

// Lots of Pism objects are returned by passing a reference to a pointer. Use 
//
// %Pism_pointer_reference_typemaps(MyType)
// 
// to declare typemaps so that an argument MyType *&OUTPUT should be treated
// as  an output (so no input variable shows up on the Python side, and
// the output shows up as an output variable).  To apply the typemap
// to all arguments of this type, use 
// %apply MyType *& OUTPUT { MyType *&}
// or use %Pism_pointer_reference_is_always_ouput(MyType) in the first place

%define %Pism_pointer_reference_typemaps(TYPE)
%typemap(in, numinputs=0,noblock=1) TYPE *& OUTPUT (TYPE *temp) {
    $1 = &temp;
}
%typemap(argout,noblock=1) TYPE *& OUTPUT
{
    %append_output(SWIG_NewPointerObj(%as_voidptr(*$1), $*descriptor, 0 | %newpointer_flags));
};
%enddef

%define %Pism_pointer_reference_is_always_output(TYPE)
%Pism_pointer_reference_typemaps(TYPE);
%apply TYPE *& OUTPUT { TYPE *&}
%enddef

/* Type maps for treating pointer-to-pointer arguments as output. */
%define %Pism_pointer_pointer_typemaps(TYPE)
%typemap(in, numinputs=0,noblock=1) TYPE ** OUTPUT (TYPE *temp) {
    $1 = &temp;
}
%typemap(argout,noblock=1) TYPE ** OUTPUT
{
    %append_output(SWIG_NewPointerObj(%as_voidptr(*$1), $*descriptor, 0 | %newpointer_flags));
};
%enddef

/* Tell SWIG that all pointer-to-pointer arguments are output. */
%define %Pism_pointer_pointer_is_always_output(TYPE)
%Pism_pointer_pointer_typemaps(TYPE);
%apply TYPE ** OUTPUT { TYPE **}
%enddef

/* Type map for treating reference arguments as output. */
%define %Pism_reference_output_typemaps(TYPE)
%typemap(in, numinputs=0,noblock=1) TYPE & OUTPUT (TYPE temp) {
    $1 = &temp;
}
%typemap(argout,noblock=1) TYPE & OUTPUT
{
    %append_output(SWIG_NewPointerObj(%as_voidptr($1), $descriptor, 0 | %newpointer_flags));
};
%enddef

/* Tell SWIG that reference arguments are always output. */
%define %Pism_reference_is_always_output(TYPE)
%Pism_reference_output_typemaps(TYPE);
%apply TYPE & OUTPUT { TYPE &}
%enddef

%typemap(in, numinputs=0,noblock=1) bool & OUTPUT (bool temp = false) {
    $1 = &temp;
}

%typemap(argout,noblock=1) bool & OUTPUT
{
    %append_output(SWIG_From(bool)(*$1));
};

%typemap(in, numinputs=0,noblock=1) PETScInt & OUTPUT (PETScInt temp) {
    $1 = &temp;
}

%typemap(argout,noblock=1) PETScInt & OUTPUT
{
    %append_output(SWIG_From(int)(*$1));
};

%typemap(in, numinputs=0,noblock=1) std::string& result (std::string temp) {
    $1 = &temp;
}

%typemap(in, numinputs=0,noblock=1) std::string& OUTPUT (std::string temp) {
    $1 = &temp;
}

%typemap(argout,noblock=1) std::string & OUTPUT
{
    %append_output(SWIG_FromCharPtr((*$1).c_str()));
}

%apply std::string &OUTPUT { std::string &result}

%typemap(in, numinputs=0,noblock=1) std::vector<int> & OUTPUT (std::vector<int> temp) {
    $1 = &temp;
}

%typemap(argout,noblock=1) std::vector<int> & OUTPUT
{
    int len;
    len = $1->size();
    $result = PyList_New(len);
     int i;
     for(i=0; i<len; i++)
     {
         PyList_SetItem($result, i, PyInt_FromLong((*$1)[i]));
     }
}

%typemap(in, numinputs=0,noblock=1) std::vector<double> & OUTPUT (std::vector<double> temp) {
    $1 = &temp;
}

%typemap(argout,noblock=1) std::vector<double> & OUTPUT
{
    int len;
    len = $1->size();
    $result = PyList_New(len);
     int i;
     for(i=0; i<len; i++)
     {
         PyList_SetItem($result, i, PyFloat_FromDouble((*$1)[i]));
     }
}

%apply std::vector<int> & OUTPUT {std::vector<int> &result};
%apply std::vector<double> & OUTPUT {std::vector<double> &result};
%apply std::vector<std::string> & OUTPUT {std::vector<std::string> & result};
 
%apply int &OUTPUT {int &result};
%apply int *OUTPUT {int *out_mask};

%apply double & OUTPUT {double & result};
%apply double & OUTPUT {double & out};
%apply double * OUTPUT {double * result};
%apply bool & OUTPUT {bool & is_set, bool & result, bool & flag, bool & success};

%Pism_pointer_pointer_is_always_output(pism::IceFlowLaw)

%include pism_options.i

// The varargs to verbPrintf aren't making it through from python.  But that's ok: we'd like
// to extend the printf features of verbPrintf to include python's formatting for objects.
// So we rename verbPrintf here and call it (without any varargs) from a python verbPrintf.
%rename(_verbPrintf) verbPrintf;

// FIXME: the the following code blocks there are explicit calls to Py????_Check.  There seems to 
// be a more elegant solution using SWIG_From(int) and so forth that I'm not familiar with.  The
// following works for now.

// The SWIG built-in typecheck for a const char [] (used, e.g., with overloaded methods) checks that 
// the string is zero length. So we have this bug fix from SWIG developer William Fulton here.
%typemap(typecheck,noblock=1,precedence=SWIG_TYPECHECK_STRING, fragment="SWIG_AsCharPtrAndSize") const char[] {
 int res = SWIG_AsCharPtrAndSize($input, 0, NULL, 0);
 $1 = SWIG_CheckState(res);
}


// Tell SWIG that the following variables are truly constant
%immutable pism::PISM_Revision;
%immutable pism::PISM_DefaultConfigFile;

/* PISM header with no dependence on other PISM headers. */
%include "enthalpyConverter.hh"
%ignore pism::Vector2::operator=;
%include "Vector2.hh"

%ignore pism::Unit::operator=;
%feature("valuewrapper") pism::UnitSystem;
%feature("valuewrapper") pism::Unit;

%include "PISMUnits.hh"
%include pism_DM.i
/* End of independent PISM classes. */

%include pism_PIO.i

/* make sure PIO.i is included before NCVariable.hh */
%include pism_NCVariable.i
%include "PISMConfig.hh"
%include "pism_const.hh"

%include pism_IceModelVec.i

%include pism_Vars.i

%include pism_Timeseries.i

%include pism_IceGrid.i

%include "PISMDiagnostic.hh"
%include "PISMComponent.hh"
%include "basal_resistance.hh"
%include "rheology/flowlaws.hh"

%include "flowlaw_factory.hh"

%include pism_ColumnSystem.i

/* SSAForwardRunFromInputFile sets up a yield stress model, which
 * requires a hydrology model.
 */
%{
using pism::StressBalance;
%}
%include "PISMHydrology.hh"

%include "Mask.hh"
%include "pism_python.hh"
%include "PISMYieldStress.hh"
%include "PISMMohrCoulombYieldStress.hh"
%include "PISMTime.hh"

%include pism_SSA.i

%include pism_SIA.i

/* The regional model implements some classes derived from SSAFD and
 * SIAFD, so this %include has to appear after %including the rest of
 * PISM's stress balance headers.
 */
%{
#include "regional/regional.hh"
%}
%include "regional/regional.hh"

%include pism_inverse.i


