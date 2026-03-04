#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include <tpp_interface.hpp> 

namespace py = pybind11;
using namespace tpp;

// Helper function to convert dpoint<double, 2> to Python list
py::list point_to_list(const reviver::dpoint<double, 2>& point) {
    py::list lst;
    lst.append(point[0]);
    lst.append(point[1]);
    return lst;
}

// Helper function to convert vector<dpoint<double, 2>> to Python list
py::list points_to_list(const std::vector<reviver::dpoint<double, 2>>& points) {
    py::list lst;
    for(auto& point : points) {
        lst.append(point_to_list(point));
    }
    return lst;
}

// Helper function to convert dpoint<double, 4> to Python list
py::list point4_to_list(const reviver::dpoint<double, 4>& point) {
    return py::make_tuple(point[0], point[1], point[2], point[3]);
}

// Helper function to convert Python list/tuple to dpoint<double, 2>
reviver::dpoint<double, 2> list_to_point(const py::list& lst) {
    if (lst.size() != 2) throw py::value_error("Point must have 2 coordinates");
    return reviver::dpoint<double, 2>(lst[0].cast<double>(), lst[1].cast<double>());
}

// Helper function to convert Python list/tuple to vector<dpoint<double, 2>>
std::vector<reviver::dpoint<double, 2>> list_to_points(const py::list& lst) {
    std::vector<reviver::dpoint<double, 2>> vec;
    for(auto& elem: lst) {
        vec.push_back(list_to_point(elem.cast<py::list>()));
    }
    return vec;
}

// Helper function to convert Python list/tuple to dpoint<double, 4>
reviver::dpoint<double, 4> list_to_point4(const py::list& lst) {
    if (lst.size() != 4) throw py::value_error("Point4 must have 4 coordinates");
    double arr[4] = { lst[0].cast<double>(), lst[1].cast<double>(),
                      lst[2].cast<double>(), lst[3].cast<double>() };
    return reviver::dpoint<double, 4>(arr);
}

// Helper class for Voronoi vertex iteration
//   OPEN TODO:: move to Triangle++ !!!
struct TRPP_LIB_EXPORT VoronoiVertexList
{
    VoronoiVertexList(Delaunay* triangulator) : m_delaunay(triangulator) {}

    struct VertexListIterator : public VoronoiVertexIterator
    {
        VertexListIterator(VoronoiVertexIterator vit) : VoronoiVertexIterator(vit) {}

        VertexListIterator operator++() {
        return VoronoiVertexIterator::operator++();
        }

        const VoronoiVertexIterator& operator*() const {
        return *this;
        }         
    };

    VertexListIterator begin() { return m_delaunay->vvbegin(); }
    VertexListIterator end() { return m_delaunay->vvend(); }

private:
    Delaunay* m_delaunay;
};

// Helper class for Voronoi edge iteration
//   OPEN TODO:: move to Triangle++ !!!
struct TRPP_LIB_EXPORT VoronoiEdgeList
{
    VoronoiEdgeList(Delaunay* triangulator) : m_delaunay(triangulator) {}

    struct VertexListIterator : public VoronoiEdgeIterator
    {
        VertexListIterator(VoronoiEdgeIterator vit) : VoronoiEdgeIterator(vit) {}

        VertexListIterator operator++() {
            return VoronoiEdgeIterator::operator++();
        }

        const VoronoiEdgeIterator& operator*() const {
            return *this;
        }         
    };

    VertexListIterator begin() { return m_delaunay->vebegin(); }
    VertexListIterator end() { return m_delaunay->veend(); }

private:
    Delaunay* m_delaunay;
};


// Helper for iterators - move to Triangle++ !!!!
//  --> not working! 
//  OPEN TODO::: modernize ietartors in Triangle++!
namespace std {
template<>
struct iterator_traits<FaceIterator> {
    using difference_type   = std::ptrdiff_t;
    using value_type        = FaceIterator::Face;
    using pointer           = FaceIterator::Face*;
    using reference         = FaceIterator::Face&;
    using iterator_category = std::random_access_iterator_tag;
};
}


PYBIND11_MODULE(triangle_ppy, m) {

    // Bind enums
    py::enum_<DebugOutputLevel>(m, "DebugOutputLevel")
        .value("Nothing", DebugOutputLevel::None)
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
        .def(py::init([](const py::list& lst, bool enableMeshIdx) { 
             return std::make_unique<Delaunay>(list_to_points(lst), enableMeshIdx);
            }),
            py::arg("points") = py::list(),
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
             py::arg("min_angle"), py::arg("max_area"))             
        .def("set_min_angle", &Delaunay::setMinAngle, py::arg("angle"))
        .def("set_max_area", &Delaunay::setMaxArea, py::arg("area"))
        .def("remove_quality_constraints", &Delaunay::removeQualityConstraints)
        
        .def("set_segment_constraint",
             [](Delaunay& d, const py::list& segments) { 
                d.setSegmentConstraint(list_to_points(segments));
             },
            py::arg("segments"))
        .def("set_segment_constraint",
             py::overload_cast<const std::vector<int>&, DebugOutputLevel>(&Delaunay::setSegmentConstraint),
             py::arg("segment_point_indexes"), py::arg("trace_level") = DebugOutputLevel::None)

        .def("use_convex_hull_with_segments", &Delaunay::useConvexHullWithSegments, py::arg("use_convex_hull"))

        .def("set_holes_constraint", 
             [](Delaunay& d, const py::list& holes) { 
                d.setHolesConstraint(list_to_points(holes));
             },            
            py::arg("holes"))

        .def("set_regions_constraint",
             [](Delaunay& d, const py::list& regions, const py::list& areas) { 
                std::vector<float> fvec;
                for(auto& a: areas) {
                    fvec.push_back(a.cast<float>());
                }
                d.setRegionsConstraint(list_to_points(regions), fvec);
             },            
             py::arg("regions"), py::arg("areas"))

        .def("check_constraints", 
               [](Delaunay& d) { 
                bool possible;
                return d.checkConstraints(possible);
             })
        .def("check_constraints_relaxed", 
               [](Delaunay& d) { 
                bool relaxed = true;
                return d.checkConstraintsOpt(relaxed);
             })
        .def("get_min_angle_boundaries",
               [](Delaunay& d) -> py::tuple { 
                float guaranteed;
                float possible;
                d.getMinAngleBoundaries(guaranteed, possible);
                return py::make_tuple(guaranteed, possible);
             })                    

        // Results API
        .def("has_triangulation", &Delaunay::hasTriangulation)
        .def("edge_count", &Delaunay::edgeCount)
        .def("triangle_count", &Delaunay::triangleCount)
        .def("vertice_count", &Delaunay::verticeCount)
        .def("hull_size", &Delaunay::hullSize)
        .def("hole_count", &Delaunay::holeCount)

        .def("get_min_max_points", 
               [](Delaunay& d) -> py::list { 
                double minX, minY, maxX, maxY;
                d.getMinMaxPoints(minX, minY, maxX, maxY);
                auto minPt = point_to_list(Delaunay::Point(minX, minY));
                auto maxPt = point_to_list(Delaunay::Point(maxX, maxY));
                return py::make_tuple(minPt, maxPt);
             })
             
        .def("faces", &Delaunay::faces) // Assumes FacesList has its own binding
        .def("vertices", &Delaunay::vertices) // Assumes VertexList has its own binding
        
        .def("voronoi_point_count", &Delaunay::voronoiPointCount)
        .def("voronoi_edge_count", &Delaunay::voronoiEdgeCount)
        
        .def("voronoi_vertices",
                [](Delaunay& d) { 
                    return VoronoiVertexList(&d);
             })
        .def("voronoi_edges",
                [](Delaunay& d) { 
                    return VoronoiEdgeList(&d);
             })


        // OPEN TODO::: ----

        .def("mesh", &Delaunay::mesh) // Assumes TriangulationMesh has its own binding
        
        
        // OPEN TODO::: ----
        
        .def("point_at_vertex_id", &Delaunay::pointAtVertexId,
             py::return_value_policy::reference)


        // File I/O API
        .def("save_points", &Delaunay::savePoints, py::arg("file_path"))
        .def("save_segments", &Delaunay::saveSegments, py::arg("file_path"))
        .def("write_off", &Delaunay::writeoff, py::arg("fname"))

        .def("read_points",
             [](Delaunay& d, const std::string& filePath, py::list& points) -> bool { 
                std::vector<reviver::dpoint<double, 2>> vec;
                if(d.readPoints(filePath, vec)) {
                    //points = points_to_list(vec); --> not working !!!                    
                    auto tmp = points_to_list(vec);
                    for(auto& point: tmp)
                        points.append(point);                    
                    return true;
                }
                return false;
             },
             py::arg("file_path"), py::arg("points"))


        // OPEN TODO::: ----

        .def("read_segments", &Delaunay::readSegments,
             py::arg("file_path"), py::arg("points"), py::arg("segment_endpoints"),
             py::arg("hole_markers"), py::arg("region_constr"),
             py::arg("duplicate_point_count") = nullptr,
             py::arg("trace_level") = DebugOutputLevel::None)
             
        .def("enable_file_io_trace", &Delaunay::enableFileIOTrace, py::arg("enable"))
        ;


    // Bind OrderPoints struct
    py::class_<Delaunay::OrderPoints>(m, "OrderPoints")
        .def(py::init<>())
        .def("__call__", &Delaunay::OrderPoints::operator(),
             py::arg("lhs"), py::arg("rhs"))
        ;

    // Iterate over faces (i.e. oriented triangles)
    py::class_<FaceIterator>(m, "FaceIterator")
        .def(py::init<>())
        ;

    py::class_<FacesList>(m, "FacesList")
        .def(py::init<Delaunay*>())
        .def("__iter__", [](FacesList &self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())

        .def("__getitem__", [](FacesList &self, size_t index) -> FaceIterator::Face {
            auto iter = self.begin();

            //std::advance(iter, index); --> not yet working, OPEN TODO:: !!!!
            //  -- OPEN TODO::: add direct indexing in Triangle++ ????
            for(size_t i = 0; i < index; ++i)
                iter = iter++;

            if(iter == self.end())
                throw std::out_of_range("Index out of range");

            return *iter;
        })

        .def("__len__", [](FacesList &self) -> size_t {
            //return self.m_delaunay->triangleCount();
            //return static_cast<size_t>(std::distance(self.begin(), self.end()));

            // OPEN TODO::: !!!
            throw std::runtime_error("Not implemented");
            return 0;
        })

        .def("as_iterator", [](FacesList &self) -> FaceIterator {
            return self.begin();
        })        
        ;

    py::class_<FaceIterator::Face>(m, "Face")
        .def(py::init<FaceIterator*>())        
        .def("org",
             [](FaceIterator::Face& f) -> py::tuple{ 
                Delaunay::Point point;
                int idx = f.Org(&point);
                return py::make_tuple(point_to_list(point), idx);
             })
        .def("org_idx",
             [](FaceIterator::Face& f) -> int { 
                int idx = f.Org();
                return idx;
             })   
        .def("org_mesh_idx",
             [](FaceIterator::Face& f) -> int{ 
                Delaunay::Point point;
                int meshIdx;
                f.Org(point, meshIdx);
                return meshIdx;
             })                     
        .def("dest",
             [](FaceIterator::Face& f) -> py::tuple{ 
                Delaunay::Point point;
                int idx = f.Dest(&point);
                return py::make_tuple(point_to_list(point), idx);
             })             
        .def("dest_idx",
             [](FaceIterator::Face& f) -> int { 
                int idx = f.Dest();
                return idx;
             })
        .def("dest_mesh_idx",
             [](FaceIterator::Face& f) -> int{ 
                Delaunay::Point point;
                int meshIdx;
                f.Dest(point, meshIdx);
                return meshIdx;
             })                
        .def("apex",
             [](FaceIterator::Face& f) -> py::tuple{ 
                Delaunay::Point point;
                int idx = f.Apex(&point);
                return py::make_tuple(point_to_list(point), idx);
             })
        .def("apex_idx",
             [](FaceIterator::Face& f) -> int { 
                int idx = f.Apex();
                return idx;
             })
        .def("apex_mesh_idx",
             [](FaceIterator::Face& f) -> int{ 
                Delaunay::Point point;
                int meshIdx;
                f.Apex(point, meshIdx);
                return meshIdx;
             })                
        .def("area",
            py::overload_cast<>(&FaceIterator::Face::area, py::const_))

        // OPEN TODO::: ----
        .def("is_ghost",
             [](FaceIterator::Face& f) -> bool {
                // OPEN TODO::: !!!  ----

                //return f.m_iter->isGhost();
                return false; 
             })   
        ;

    // Iterate over triangle vertices        
    py::class_<VertexList>(m, "VertexList")
        .def(py::init<Delaunay*>())
        .def("__iter__", [](VertexList &self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())        
        ;

    py::class_<VertexIterator>(m, "VertexIterator")
        .def(py::init<>())
        .def("point",
             [](VertexIterator& v) -> py::list{ 
                return point_to_list(*v);
             })
        .def("x", &VertexIterator::x)
        .def("y", &VertexIterator::y)
        .def("vertex_id", &VertexIterator::vertexId)
        ;

    // Iterate over Voronoi vertices        
    py::class_<VoronoiVertexList>(m, "VoronoiVertexList")
        .def(py::init<Delaunay*>())
        .def("__iter__", [](VoronoiVertexList &self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())        
        ;

    py::class_<VoronoiVertexIterator>(m, "VoronoiVertexIterator")
        .def(py::init<>())  
        .def("point",
             [](VoronoiVertexIterator& v) -> py::list{ 
                return point_to_list(*v);
             })
        ;

    // Iterate over Voronoi edges
    py::class_<VoronoiEdgeList>(m, "VoronoiEdgeList")
        .def(py::init<Delaunay*>())
        .def("__iter__", [](VoronoiEdgeList &self) {
            return py::make_iterator(self.begin(), self.end());
        }, py::keep_alive<0, 1>())        
        ;

    py::class_<VoronoiEdgeIterator>(m, "VoronoiEdgeIterator")
        .def(py::init<>())
        .def("start_point_id",
             [](VoronoiEdgeIterator& v) -> int { 
                return v.startPointId();
             })
        .def("end_point_id",
             [](VoronoiEdgeIterator& v) -> py::tuple { 
                Delaunay::Point normvec;
                int idx = v.endPointId(normvec);
                return py::make_tuple(idx, point_to_list(normvec));
             })
        .def("org",
             [](VoronoiEdgeIterator& v) -> py::list { 
                return point_to_list(v.Org());
             })
        .def("dest",
             [](VoronoiEdgeIterator& v) -> py::tuple { 
                bool finiteEdge;
                py::list pt = point_to_list(v.Dest(finiteEdge));
                return py::make_tuple(pt, finiteEdge);
             })
        ;


    // OPEN TODO::: ----

    py::class_<TriangulationMesh>(m, "TriangulationMesh")
        .def(py::init<Delaunay*>())
        .def("opposite", 
             [](TriangulationMesh& mesh, FaceIterator& f) -> FaceIterator { 
                return mesh.Sym(f);
             }, 
             py::arg("face"))
        
            // Add methods as needed
        ;

}
