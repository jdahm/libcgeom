#include "par/environment.hpp"

namespace par
{

namespace detail
{

env::initializer::initializer() {
        int thread_mode;
        MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &thread_mode);
        // MPI_Init(0, 0);
        float ver;
        Zoltan_Initialize(0, 0, &ver);
}

env::initializer::~initializer() {
        MPI_Finalize();
}

env::env() : init(), comm_world_(MPI_COMM_WORLD), comm_self_(MPI_COMM_SELF) { }

threading_modes env::threading_mode() const {
        int provided;
        MPI_Query_thread(&provided);
        switch (provided) {
        case MPI_THREAD_SINGLE:
                return threading_modes::single;
        case MPI_THREAD_FUNNELED:
                return threading_modes::funneled;
        case MPI_THREAD_SERIALIZED:
                return threading_modes::serialized;
        case MPI_THREAD_MULTIPLE:
                return threading_modes::multiple;
        default:
                return threading_modes::single;
        }
}

bool env::is_thread_main() const {
        int res;
        MPI_Is_thread_main(&res);
        return static_cast<bool>(res);
}

const communicator & env::comm_world() const { return comm_world_; }

const communicator & env::comm_self() const { return comm_self_; }

std::string env::processor_name() const {
        char name[MPI_MAX_PROCESSOR_NAME];
        int len;
        MPI_Get_processor_name(name, &len);
        return std::string(name);
}

double env::wtime() const { return MPI_Wtime(); }

void env::buffer_attach(void *buff, int size) const {
        MPI_Buffer_attach(buff, size);
}

std::pair<void *, int> env::buffer_detach() const {
        void *buff;
        int size;
        MPI_Buffer_detach(&buff, &size);
        return std::make_pair(buff, size);
}

const env & get_env() {
        static env instance;
        return instance;
}

} // namespace detail

constexpr int any_source() {
        return MPI_ANY_SOURCE;
}

constexpr int proc_null() {
        return MPI_PROC_NULL;
}

constexpr int undefined() {
        return MPI_UNDEFINED;
}

constexpr int root() {
        return MPI_ROOT;
}

threading_modes threading_mode() {
        return detail::get_env().threading_mode();
}

bool is_thread_main() {
        return detail::get_env().is_thread_main();
}

bool is_root() {
        return detail::get_env().comm_world().rank() == 0;
}

int size() {
        return detail::get_env().comm_world().size();
}

const communicator & comm_world() {
        return detail::get_env().comm_world();
}

const communicator & comm_self() {
        return detail::get_env().comm_self();
}

std::string processor_name() {
        return detail::get_env().processor_name();
}

double wtime() {
        return detail::get_env().wtime();
}

void buffer_attach(void *buff, int size) {
        return detail::get_env().buffer_attach(buff, size);
}

std::pair<void *, int> buffer_detach() {
        return detail::get_env().buffer_detach();
}

} // namespace par
