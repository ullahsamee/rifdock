#ifndef INCLUDED_scheme_nest_storage_policy_HH
#define INCLUDED_scheme_nest_storage_policy_HH


namespace scheme {
namespace nest {


	///////////////////////////////////////////////////////////////////////////////////////////////
	/// storage policies
	///////////////////////////////////////////////////////////////////////////////////////////////

	///@brief Storage Policy Class, store by value
	///@detail this is a sortof DUMMY baseclass to define the StoragePolicy concept
	///  also a placeholder for doxygen 
	///@tparam Value type of value that is stored
	template< class Value >
	class StoragePolicy {
		///@return reference to the stored value
		Value const & value() const;
	protected:
		///@return nonconst reference to stored Value
		Value & nonconst_value();
	};

	/// @brief Storage Policy Class, store by value
	/// @tparam Value ValueType stored
	template< class Value >
	struct StoreValue : public StoragePolicy<Value> {
		///@return const reference to stored Value
		Value const & value() const { return value_; }
	protected:
		///@return nonconst reference to stored Value
		Value & nonconst_value() { return value_; }
		Value value_;
		~StoreValue(){}
	};

	/// @brief Store-by-pointer policy
	/// @tparam Value ValueType stored
	/// @note addes the ability to set the pointer
	template< class Value >
	struct StorePointer : public StoragePolicy<Value> {
		///@return const reference to stored Value
		Value const & value() const { return *value_; }
		/// @brief switch the pointer the this policy manages
		/// @param new_pointer 
		void set_pointer(Value * new_pointer) { value_ = new_pointer; }
	protected:
		///@return nonconst reference to stored Value
		Value & nonconst_value() { return *value_; }
		Value * value_;
		~StorePointer(){}
	};



}
}

#endif
