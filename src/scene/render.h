#pragma once

#include "core/lajolla.h"
#include "scene/image.h"
#include <memory>

struct Scene;

Image3 render(const Scene &scene);
