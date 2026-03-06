from triangle_ppy import Delaunay, DebugOutputLevel

points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]

d = Delaunay(points)
d.set_quality_constraints(min_angle=20.0, max_area=0.1)

d.triangulate(quality=True, trace_level=DebugOutputLevel.Info)

print(f"Triangle count: {d.triangle_count()}")

faces = d.faces() # faces are _oriented_ triangles!

for idx, f in enumerate(faces):
    point1 = f.org_pt()
    point2 = f.dest_pt()
    point3 = f.apex_pt()

    print(f"Triangle {idx}: {point1} -> {point2} -> {point3}")

try:  
    d.save_points("./min_trg_points.node")
    print(f"Resulting points saved to file") 
except:
    print(f"Error - cannot save points to file!") 
