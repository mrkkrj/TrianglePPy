#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>
#include "Delaunay.h" // Assume this includes the provided C++ class

namespace py = pybind11;
using namespace tpp;

// Helper function to convert dpoint<double, 2> to Python list
py::list point_to_list(const reviver::dpoint<double, 2>& point) {
    return py::make_tuple(point.x, point.y);
}

// Helper function to convert dpoint<double, 4> to Python list
py::list point4_to_list(const reviver::dpoint<double, 4>& point) {
    return py::make_tuple(point.x, point.y, point.z, point.w);
}

// Helper function to convert Python list/tuple to dpoint<double, 2>
reviver::dpoint<double, 2> list_to_point(const py::list& lst) {
    if (lst.size() != 2) throw py::value_error("Point must have 2 coordinates");
    return reviver::dpoint<double, 2>(lst[0].cast<double>(), lst[1].cast<double>());
}

// Helper function to convert Python list/tuple to dpoint<double, 4>
reviver::dpoint<double, 4> list_to_point4(const py::list& lst) {
    if (lst.size() != 4) throw py::value_error("Point4 must have 4 coordinates");
    return reviver::dpoint<double, 4>(
        lst[0].cast<double>(), lst[1].cast<double>(),
        lst[2].cast<double>(), lst[3].cast<double>()
    );
}

PYBIND11_MODULE(delaunay, m) {
    // Bind enums
    py::enum_<DebugOutputLevel>(m, "DebugOutputLevel")
        .value("None", DebugOutputLevel::None)
        .value("Info", DebugOutputLevel::Info)
        .value("Vertex", DebugOutputLevel::Vertex)
        .value("Debug", DebugOutputLevel::Debug)
        .export_values();

    py::enum_<AlgorithmType>(m, "AlgorithmType")
        .value("DivideConquer", AlgorithmType::DivideConquer)
        .value("Incremental", AlgorithmType::Incremental)
        .value("Sweepline", AlgorithmType::Sweepline)
        .export_values();

    // Bind the Delaunay class
    py::class_<Delaunay>(m, "Delaunay")
        // Constructor
        .def(py::init<const std::vector<reviver::dpoint<double, 2>>&, bool>(),
             py::arg("points") = std::vector<reviver::dpoint<double, 2>>(),
             py::arg("enable_mesh_indexing") = false)

        // Main API
        .def("triangulate", py::overload_cast<bool, DebugOutputLevel>(&Delaunay::Triangulate),
             py::arg("quality") = false, py::arg("trace_level") = DebugOutputLevel::None)
        .def("triangulate", py::overload_cast<DebugOutputLevel>(&Delaunay::Triangulate),
             py::arg("trace_level"))
        .def("triangulate_conf", py::overload_cast<bool, DebugOutputLevel>(&Delaunay::TriangulateConf),
             py::arg("quality") = false, py::arg("trace_level") = DebugOutputLevel::None)
        .def("triangulate_conf", py::overload_cast<DebugOutputLevel>(&Delaunay::TriangulateConf),
             py::arg("trace_level"))
        .def("tesselate", &Delaunay::Tesselate,
             py::arg("use_conforming_delaunay") = false, py::arg("trace_level") = DebugOutputLevel::None)
        .def("enable_mesh_index_generation", &Delaunay::enableMeshIndexGeneration)
        .def("set_algorithm", &Delaunay::setAlgorithm, py::arg("algorithm"))

        // Constraints API
        .def("set_quality_constraints", &Delaunay::setQualityConstraints,
             py::arg("angle"), py::arg("area"))
        .def("set_min_angle", &Delaunay::setMinAngle, py::arg("angle"))
        .def("set_max_area", &Delaunay::setMaxArea, py::arg("area"))
        .def("remove_quality_constraints", &Delaunay::removeQualityConstraints)
        .def("set_segment_constraint",
             py::overload_cast<const std::vector<reviver::dpoint<double, 2>>&>(&Delaunay::setSegmentConstraint),
             py::arg("segments"))
        .def("set_segment_constraint",
             py::overload_cast<const std::vector<int>&, DebugOutputLevel>(&Delaunay::setSegmentConstraint),
             py::arg("segment_point_indexes"), py::arg("trace_level") = DebugOutputLevel::None)
        .def("use_convex_hull_with_segments", &Delaunay::useConvexHullWithSegments, py::arg("use_convex_hull"))
        .def("set_holes_constraint", &Delaunay::setHolesConstraint, py::arg("holes"))
        .def("set_regions_constraint",
             py::overload_cast<const std::vector<reviver::dpoint<double, 2>>&, const std::vector<float>&>(
                 &Delaunay::setRegionsConstraint),
             py::arg("regions"), py::arg("areas"))
        .def("set_regions_constraint",
             py::overload_cast<const std::vector<reviver::dpoint<double, 4>>&>(&Delaunay::setRegionsConstraint),
             py::arg("region_constr"))
        .def("check_constraints", &Delaunay::checkConstraints, py::arg("possible"))
        .def("check_constraints_opt", &Delaunay::checkConstraintsOpt, py::arg("relaxed"))
        .def_static("get_min_angle_boundaries", &Delaunay::getMinAngleBoundaries,
                    py::arg("guaranteed"), py::arg("possible"))

        // Results API
        .def("has_triangulation", &Delaunay::hasTriangulation)
        .def("edge_count", &Delaunay::edgeCount)
        .def("triangle_count", &Delaunay::triangleCount)
        .def("vertice_count", &Delaunay::verticeCount)
        .def("hull_size", &Delaunay::hullSize)
        .def("hole_count", &Delaunay::holeCount)
        .def("get_min_max_points", &Delaunay::getMinMaxPoints,
             py::arg("min_x"), py::arg("min_y"), py::arg("max_x"), py::arg("max_y"))
        .def("faces", &Delaunay::faces) // Assumes FacesList has its own binding
        .def("vertices", &Delaunay::vertices) // Assumes VertexList has its own binding
        .def("voronoi_point_count", &Delaunay::voronoiPointCount)
        .def("voronoi_edge_count", &Delaunay::voronoiEdgeCount)
        .def("mesh", &Delaunay::mesh) // Assumes TriangulationMesh has its own binding
        .def("point_at_vertex_id", &Delaunay::pointAtVertexId,
             py::return_value_policy::reference)

        // File I/O API
        .def("save_points", &Delaunay::savePoints, py::arg("file_path"))
        .def("save_segments", &Delaunay::saveSegments, py::arg("file_path"))
        .def("write_off", &Delaunay::writeoff, py::arg("fname"))
        .def("read_points", &Delaunay::readPoints,
             py::arg("file_path"), py::arg("points"))
        .def("read_segments", &Delaunay::readSegments,
             py::arg("file_path"), py::arg("points"), py::arg("segment_endpoints"),
             py::arg("hole_markers"), py::arg("region_constr"),
             py::arg("duplicate_point_count") = nullptr,
             py::arg("trace_level") = DebugOutputLevel::None)
        .def("enable_file_io_trace", &Delaunay::enableFileIOTrace, py::arg("enable"))

        // Iterator methods (simplified, assuming iterator bindings)
        .def("fbegin", &Delaunay::fbegin)
        .def("fend", &Delaunay::fend)
        .def("vbegin", &Delaunay::vbegin)
        .def("vend", &Delaunay::vend)
        .def("vvbegin", &Delaunay::vvbegin)
        .def("vvend", &Delaunay::vvend)
        .def("vebegin", &Delaunay::vebegin)
        .def("veend", &Delaunay::veend);

    // Bind OrderPoints struct
    py::class_<Delaunay::OrderPoints>(m, "OrderPoints")
        .def(py::init<>())
        .def("__call__", &Delaunay::OrderPoints::operator(),
             py::arg("lhs"), py::arg("rhs"));

    // Placeholder bindings for iterators and other classes (to be expanded)
    py::class_<FaceIterator>(m, "FaceIterator")
        .def(py::init<>())
        // Add iterator protocol
        .def("__iter__", [](FaceIterator &it) -> FaceIterator& { return it; })
        .def("__next__", [](FaceIterator &it) {
            // Implement next logic or throw StopIteration
            throw py::stop_iteration();
        });

    py::class_<VertexIterator>(m, "VertexIterator")
        .def(py::init<>())
        .def("__iter__", [](VertexIterator &it) -> VertexIterator& { return it; })
        .def("__next__", [](VertexIterator &it) {
            throw py::stop_iteration();
        });

    py::class_<VoronoiVertexIterator>(m, "VoronoiVertexIterator")
        .def(py::init<>())
        .def("__iter__", [](VoronoiVertexIterator &it) -> VoronoiVertexIterator& { return it; })
        .def("__next__", [](VoronoiVertexIterator &it) {
            throw py::stop_iteration();
        });

    py::class_<VoronoiEdgeIterator>(m, "VoronoiEdgeIterator")
        .def(py::init<>())
        .def("__iter__", [](VoronoiEdgeIterator &it) -> VoronoiEdgeIterator& { return it; })
        .def("__next__", [](VoronoiEdgeIterator &it) {
            throw py::stop_iteration();
        });

    py::class_<TriangulationMesh>(m, "TriangulationMesh")
        .def(py::init<>())
        // Add methods as needed
        ;

    py::class_<FacesList>(m, "FacesList")
        .def(py::init<>())
        // Add methods as needed
        ;

    py::class_<VertexList>(m, "VertexList")
        .def(py::init<>())
        // Add methods as needed
        ;
}
