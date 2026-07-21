#pragma once

#include "core/lajolla.h"
#include "core/math/matrix.h"
#include "geometry/shape.h"

/// Load Mitsuba's serialized file format.
TriangleMesh load_serialized(const fs::path &filename,
                             int shape_index,
                             const Matrix4x4 &to_world);
