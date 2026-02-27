from triangle_ppy import Delaunay, DebugOutputLevel

# Create a Delaunay object with points
points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]

d = Delaunay(points, enable_mesh_indexing=True)

# Set constraints
d.set_quality_constraints(angle=20.0, area=0.1)

#d.set_segment_constraint([[0.0, 0.0], [1.0, 0.0]])
# same as this:
#d.set_segment_constraint(segment_point_indexes = [0, 1])

# Perform triangulation
d.triangulate(quality=True, trace_level=DebugOutputLevel.Info)

# Get results
print(f"Triangle count: {d.triangle_count()}")
print(f"Edge count: {d.edge_count()}")

# Save to file
d.save_points("triangulation_result.node")

# Read from file
points_from_file = []
d.read_points("triangulation_result.node", points_from_file)
#print(f"Points from file count: {len(points_from_file)}") # not working!! why???
print(f"Points from file: {points_from_file}") 
