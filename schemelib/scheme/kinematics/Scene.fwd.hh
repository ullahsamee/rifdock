#ifndef INCLUDED_kinematics_Scene_fwd_HH
#define INCLUDED_kinematics_Scene_fwd_HH

#include "scheme/types.hh"
#include <vector>

namespace scheme {
namespace kinematics {


template<
	class _Position,
	class _Index = uint64_t
>
struct SceneBase;


template<
	class _Conformation,
	class _Position,
	class _Index = uint64_t
>
struct Scene;


namespace impl {

	template< class ActorContainer, class _CacheData=uint64_t >
	struct Conformation;


}
}
}

#endif
