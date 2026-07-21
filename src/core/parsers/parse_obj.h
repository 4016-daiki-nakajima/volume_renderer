#pragma once

#include "core/lajolla.h"
#include "core/math/matrix.h"
#include "geometry/shape.h"

/// Parse Wavefront obj files. Currently only supports triangles and quads.
/// Throw errors if encountered general polygons.
TriangleMesh parse_obj(const fs::path &filename, const Matrix4x4 &to_world);
