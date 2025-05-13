from delaunay import Delaunay, DebugOutputLevel, AlgorithmType

# Create a Delaunay object with points
points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]
d = Delaunay(points, enable_mesh_indexing=True)

# Set constraints
d.set_quality_constraints(angle=20.0, area=0.1)
d.set_segment_constraint([[0.0, 0.0], [1.0, 0.0]])

# Perform triangulation
d.triangulate(quality=True, trace_level=DebugOutputLevel.Info)

# Get results
print(f"Triangle count: {d.triangle_count()}")
print(f"Vertices: {d.vertices()}")

# Save to file
d.save_points("output.node")
