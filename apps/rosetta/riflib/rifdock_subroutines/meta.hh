
    #include <riflib/rifdock_subroutines/util.hh>
    
    #include <riflib/rifdock_subroutines/hsearch_original.hh>

    #include <riflib/rifdock_subroutines/hack_pack.hh>
    #include <riflib/rifdock_subroutines/rosetta_rescore.hh>
    #include <riflib/rifdock_subroutines/compile_and_filter_results.hh>
    #include <riflib/rifdock_subroutines/output_results.hh>

#include <scheme/util/meta/util.hh>
#include <scheme/kinematics/Director.meta.hh>
#include <riflib/rifdock_typedefs.hh>
    

using ::scheme::make_shared;
using ::scheme::shared_ptr;

template <class __Director>
using _DirectorBase = shared_ptr< ::scheme::kinematics::Director<_DirectorPosition<__Director>,
    _DirectorBigIndex<__Director>,
    _DirectorIndex<__Director>
     >>;


template <class __Director>
using _SearchPointWithRots = tmplSearchPointWithRots<_DirectorBigIndex<__Director>>;

template <class __Director>
using _RifDockResult = tmplRifDockResult<_DirectorBigIndex<__Director>>;

template <class __Director>
using _SearchPoint = tmplSearchPoint<_DirectorBigIndex<__Director>>;
