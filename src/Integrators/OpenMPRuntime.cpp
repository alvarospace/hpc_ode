#include <memory>
#include <unordered_map>
#include <string>
#include <sstream>
#include <algorithm>

#include <omp.h>

#include "yaml-cpp/yaml.h"

#include "ODEIntegrator/Logger/Logger.hpp"

#include "Integrators/OpenMPRuntime.hpp"

void setOMPRuntime(YAML::Node ompConfig, std::shared_ptr<BaseLogger> logger) {
    // OMP configuration
    std::string schedule = ompConfig["schedule"]["type"].as<std::string>();
    int const scheduleChunk = ompConfig["schedule"]["chunk"].as<int>();
    int const numCpus = ompConfig["cpus"].as<int>();

    std::unordered_map<std::string, omp_sched_t> scheduleMap = {
        {"static", omp_sched_t::omp_sched_static},
        {"dynamic", omp_sched_t::omp_sched_dynamic},
        {"guided", omp_sched_t::omp_sched_guided},
        {"auto", omp_sched_t::omp_sched_auto}
    };

    std::stringstream ss;
    
    // Set numCpus or max 
    int const maxCpus = omp_get_max_threads();
    ss << "Max OpenMP threads available: " << maxCpus;
    logger->info(ss.str());
    ss.str(std::string());

    ss << "Num OpenMP threads wanted: " << numCpus;
    ss.str(std::string());
    int realCpus = std::min<int>(numCpus, maxCpus);
    omp_set_num_threads(realCpus);

    #pragma omp parallel
    {
        int const master_id = 0;
        int thread_id = omp_get_thread_num();
        if (thread_id == master_id) {
            int runningCpus = omp_get_num_threads();
            ss << "Num OpenMP threads established: " << runningCpus;
            logger->info(ss.str());
            ss.str(std::string());
        }
    }

    omp_set_schedule(scheduleMap[schedule], scheduleChunk);
    omp_sched_t realSchedule;
    int realChunk;
    omp_get_schedule(&realSchedule, &realChunk);

    std::string realScheduleStr;
    for (auto const& [key, not_needed] : scheduleMap) {
        if (scheduleMap[key] == realSchedule) {
            realScheduleStr = key;
            break;
        }
    }

    ss << "OpenMP schedule established: " << realScheduleStr;
    logger->info(ss.str());
    ss.str(std::string());

    ss << "OpenMP schedule chunk established: ";
    if (realChunk < 1) {
        ss << "auto";
    } else {
        ss << realChunk;
    }
    logger->info(ss.str());
    ss.str(std::string());
}
