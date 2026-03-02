from triangle_ppy import Delaunay, DebugOutputLevel

output_level = DebugOutputLevel.Info  #.Nothing #.Debug #.Info 

# Create a Delaunay object with points
points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]

d = Delaunay(points, enable_mesh_indexing=True)


# Set constraints
d.set_quality_constraints(min_angle=20.0, max_area=0.1)

d.set_segment_constraint([[0.0, 0.1], [1.0, 1.0]])
# same as this:
#d.set_segment_constraint(segment_point_indexes = [2, 3])

# Perform triangulation
d.triangulate(quality=True, trace_level=output_level)

# Get results
triangle_count = d.triangle_count()
edge_count = d.edge_count()
vertice_count = d.vertice_count()

if output_level != DebugOutputLevel.Nothing:
    print(f"Triangle count: {triangle_count}")
    print(f"Edge count: {edge_count}")
    print(f"Vertice count: {vertice_count}")
    # etc ...

print(f" >> Triangulation OK ... \n")    


# Iterate over triangles
faces = d.faces()
for f in faces:
    point1, input_points_idx1 = f.org()
    point2, input_points_idx2 = f.dest()
    point3, input_points_idx3 = f.apex()

    # you can ignore the index!
    pointApx = f.apex()
    pointOrg, _ = f.org()

    # or get just the index:
    orig_index = f.org_idx()

    # also the index in the resulting mesh (in iteration order!)
    orig_mesh_index = f.org_mesh_idx()

    # access data
    #  - if index == -1: an added point (i.e. a Steiner pt!)
    p1_x = points[input_points_idx1][0] if input_points_idx1 != -1 else point1[0]
    p1_y = points[input_points_idx1][1] if input_points_idx1 != -1 else point1[1]

    p2_x = points[input_points_idx2][0] if input_points_idx2 != -1 else point2[0]
    p2_y = points[input_points_idx2][1] if input_points_idx2 != -1 else point2[1]

    p3_x = points[input_points_idx3][0] if input_points_idx3 != -1 else point3[0]
    p3_y = points[input_points_idx3][1] if input_points_idx3 != -1 else point3[1]

    if output_level != DebugOutputLevel.Nothing:
        print(f"Triangle by pts: {point1}/{point2}/{point3} - area: {f.area()}")  
        print(f"Triangle by idx: {[p1_x, p1_y]}/{[p2_x, p2_y]}/{[p3_x, p3_y]}")
        print(f" -- Origin: index in input={orig_index}, index in mesh={orig_mesh_index}")            
    
# Iterate over vertices
vertices = d.vertices()
for v in vertices:
    pt = v.point()
    p_x = v.x()
    p_y = v.y()
    v_id = v.vertex_id()

    if output_level != DebugOutputLevel.Nothing:
        print(f"vertexID: {v_id} - at {pt}") 

print(f" >> Iteration OK ... \n")


# Save to file
d.save_points("./triangulation_result.node")

# Read from file
points = []
if(d.read_points("./triangulation_result.node", points)):    
    if output_level != DebugOutputLevel.Nothing:
        print(f"Points from file count: {len(points)}")
        print(f"Points from file: {points}") 
else:
    print(f"Cannot read points from file: {points}") 

print(f" >> File I/O OK ... \n")


# Perform tesselation
d.tesselate(use_conforming_delaunay=False, trace_level=output_level)

# Get results
point_count = d.voronoi_point_count()
edge_count = d.voronoi_edge_count()

if output_level != DebugOutputLevel.Nothing:
    print(f"Voronoi point count: {point_count}")
    print(f"Voronoi edge count: {edge_count}")

print(f" >> Tesselation OK ... \n")


# Iterate over the diagram
vertices = d.voronoi_vertices()
for v in vertices:
    pt = v.point()

    if output_level != DebugOutputLevel.Nothing:
        print(f"Voronoi point: {pt}") 

edges = d.voronoi_edges()
for e in edges:
    orig = e.org()
    dest, is_finite_edge = e.dest()

    if output_level != DebugOutputLevel.Nothing:
        print(f"Voronoi edge: {orig} -> {dest if is_finite_edge else ['inf', 'inf']}") 

print(f" >> Voronoi Iteration OK ... \n")


# --- more ????


print(f"Finished ---> \n")