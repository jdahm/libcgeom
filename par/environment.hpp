#ifndef PAR_ENVIRONMENT_HPP
#define PAR_ENVIRONMENT_HPP

// Based on mpl bindings: https://github.com/rabauke/mpl

#include <mpi.h>
#include <string>
#include <utility>

#include "par/traits.hpp"
#include "par/communicator.hpp"

namespace par
{
namespace detail
{

class env {
	class initializer {
	public:
                initializer();
                ~initializer();
	};

        initializer init;
        communicator comm_world_, comm_self_;
public:
        env();

	env(const env&) = delete;
	env& operator=(const env&) = delete;

	threading_modes threading_mode() const;

	bool is_thread_main() const;

	const communicator & comm_world() const;
	const communicator & comm_self() const;

	std::string processor_name() const;

        double wtime() const;

	void buffer_attach(void*, int) const;
	std::pair<void*, int> buffer_detach() const;
};

const env& get_env();

} // namespace detail

constexpr int any_source();
constexpr int proc_null();
constexpr int undefined();
constexpr int root();

threading_modes threading_mode();

bool is_thread_main();
bool is_root();

int size();

const communicator & comm_world();
const communicator & comm_self();

std::string processor_name();

double wtime();

void buffer_attach(void *buff, int size);
std::pair<void *, int> buffer_detach();

} // namespace par

#endif
