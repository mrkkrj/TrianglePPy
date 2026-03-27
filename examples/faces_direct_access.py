from triangle_ppy import Delaunay, DebugOutputLevel, enable_debug_trace

output_level = DebugOutputLevel.Info  #.Nothing #.Debug #.Info 

# Create a Delaunay object with points
points = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [1.0, 1.0]]
d = Delaunay(points, enable_mesh_indexing=True)

# Set constraints
d.set_quality_constraints(min_angle=20.0, max_area=0.1)
d.set_segment_constraint([[0.0, 0.1], [1.0, 1.0]])

ok = d.check_constraints()
assert ok, "Constraints should be OK"

print(f" >> Constraints set OK ... \n")  

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
for i, f in enumerate(faces):
    face_area = f.area()
    if output_level != DebugOutputLevel.Nothing:
        print(f"Triangle {i} - area: {face_area}")  

print(f" >> Face iteration OK ... \n")    


#####################################
#
#   ---> WIP !!!!
#
######################################

enable_debug_trace(True)

# faces can be also adressed directly
#  -- caution: horrible performance for the moment!!! WIP...
faces_ct = len(faces)
if output_level != DebugOutputLevel.Nothing:
    print(f"Face count: {faces_ct}")
      
assert faces_ct == triangle_count # as faces are *oriented* triangles!

print(f"Face direct access TEST:: face[0] ...") 
first_face = faces[0]
print(f"Face direct access TEST:: face[1] ...") 
second_face = faces[1]
print(f"Face direct access TEST:: face[15] ...") 
second_face = faces[15]

# OPEN TODO:::
# -- not yet working, GC removes face object!!!

print(f"Face direct access TEST:: use face[0] ...") 

first_face_area = faces[0].area() # crashes ??!!!
print(f" --- first_face_area: {first_face_area}")  

print(f"Face direct access TEST:: use (saved) first_face=face[0] ...") 

first_face_area = first_face.area()
print(f" --- first_face_area: {first_face_area}")  


if faces_ct > 0:
    first_face = faces[0]
    if output_level != DebugOutputLevel.Nothing:
        # not working!
        #print(f"First face: {first_face.org_pt()}/{first_face.dest_pt()}/{first_face.apex_pt()} - area: {first_face.area()}")  
        print(f"First face: {first_face}")  
        #print(f"First face - area: {first_face.area()}")  # crashes!!! ?????

print(f"Faces[] OK: {faces_ct}")  

for ff in faces:
    print(f"face: {ff} area {ff.area()}")  

print(f"Faces area loop OK: {faces_ct}")  

print(f" Trying to use first face of {faces_ct} ........")  
print(f"First face - area: {first_face.area()}")  # crashes!!! ?????
print(f"First face OK: {faces_ct} \n ---- \n")  

# ....


print(f"Finished ---> \n")