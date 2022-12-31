#pragma once

#include "yaml-cpp/yaml.h"

#include "ODEIntegrator/Logger/BaseLogger.hpp"

void setOMPRuntime(YAML::Node ompConfig, std::shared_ptr<BaseLogger> logger);
