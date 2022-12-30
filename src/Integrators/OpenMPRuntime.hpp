#pragma once

#include <unordered_map>
#include <string>

#include <omp.h>

#include "yaml-cpp/yaml.h"

// TODO: finish this runtime omp configuration
void setOMPRuntime(YAML::Node ompConfig) {
    // OMP configuration
    std::string const schedule = ompConfig["schedule"]["type"].as<std::string>();
    int const scheduleChunk = ompConfig["schedule"]["chunk"].as<int>();
    int const numCpus = ompConfig["cpus"].as<int>();

    std::unordered_map<std::string, omp_sched_t> const scheduleMap = {
        {"static", omp_sched_t::omp_sched_static},
        {"dynamic", omp_sched_t::omp_sched_dynamic},
        {"guided", omp_sched_t::omp_sched_guided},
        {"auto", omp_sched_t::omp_sched_auto}
    };


}