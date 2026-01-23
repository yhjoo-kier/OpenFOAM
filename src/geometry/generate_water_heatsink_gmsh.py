#!/usr/bin/env python3
import gmsh
import sys
import math

def create_geometry(output_file):
    gmsh.initialize()
    gmsh.model.add("heatsink_water")

    # Dimensions
    L_dom = 0.080  # 80mm
    W_dom = 0.083  # 83mm
    H_base = 0.002 # 2mm
    H_fin = 0.015  # 15mm
    
    # Fins
    D_fin = 0.006  # 6mm
    R_fin = D_fin / 2
    N_rows = 8
    N_cols = 8
    pitch = 0.010 # 10mm
    
    # Start positions (center of first fin)
    x0_fin = 0.005
    y0_fin = 0.0065
    
    # Tubes
    D_tube = 0.005 # 5mm
    R_tube = D_tube / 2
    L_tube = 0.025 # 25mm
    z_tube = 0.0075 # 7.5mm (from bottom of base)
    
    # Inlet tube
    y_inlet_center = 0.005
    # Outlet tube
    y_outlet_center = 0.078
    
    lc_solid = 0.002
    lc_fluid = 0.0015

    factory = gmsh.model.occ

    # 1. Create Base Plate (Solid)
    # Box: x, y, z, dx, dy, dz
    base_tag = factory.addBox(0, 0, 0, L_dom, W_dom, H_base)
    
    # 2. Create Fins (Solid)
    fin_tags = []
    for i in range(N_rows): # x direction
        for j in range(N_cols): # y direction
            xc = x0_fin + i * pitch
            yc = y0_fin + j * pitch
            # Cylinder: x, y, z, dx, dy, dz, r
            ft = factory.addCylinder(xc, yc, H_base, 0, 0, H_fin, R_fin)
            fin_tags.append(ft)
            
    # Fuse base and fins to create single solid domain
    solid_tag, _ = factory.fuse([(3, base_tag)], [(3, t) for t in fin_tags])
    solid_vol_tag = solid_tag[0][1]

    # 3. Create Fluid Domain
    # Main channel
    fluid_block = factory.addBox(0, 0, H_base, L_dom, W_dom, H_fin)
    
    # Inlet Tube
    inlet_tube = factory.addCylinder(0, y_inlet_center, z_tube, -L_tube, 0, 0, R_tube)
    
    # Outlet Tube
    outlet_tube = factory.addCylinder(L_dom, y_outlet_center, z_tube, L_tube, 0, 0, R_tube)
    
    # Fuse fluid parts
    fluid_parts, _ = factory.fuse([(3, fluid_block)], [(3, inlet_tube), (3, outlet_tube)])
    fluid_vol_pre = fluid_parts[0][1]

    # 4. Interface Handling
    # Use BooleanFragments to ensure conformal mesh at the interface
    ov, ov_map = factory.fragment([(3, fluid_vol_pre)], [(3, solid_vol_tag)])
    
    factory.synchronize()
    
    # Identify volumes
    # ov contains all resulting volumes [dim, tag]
    # ov_map maps input to output. 
    # Input 0 (Fluid) -> [Fluid_Part, Solid_Part]
    # Input 1 (Solid) -> [Solid_Part]
    
    # The Solid_Part is the common volume.
    solid_final_tag = ov_map[1][0][1]
    
    # The Fluid_Part is the one in Fluid map but not in Solid map
    fluid_final_tag = -1
    for t in ov_map[0]:
        if t[1] != solid_final_tag:
            fluid_final_tag = t[1]
            break
            
    if fluid_final_tag == -1:
        print("Error: Could not identify fluid volume")
        sys.exit(1)

    # 5. Physical Groups
    gmsh.model.addPhysicalGroup(3, [solid_final_tag], name="solid")
    gmsh.model.addPhysicalGroup(3, [fluid_final_tag], name="fluid")
    
    # 6. Physical Surfaces
    surfaces = gmsh.model.getEntities(dim=2)
    
    inlet_surfs = []
    outlet_surfs = []
    heat_source_surfs = []
    fluid_walls_surfs = []
    solid_walls_surfs = []
    interface_surfs = [] # Not strictly needed if OpenFOAM splits automatically, but good to check
    
    eps = 1e-4
    
    for surf in surfaces:
        tag = surf[1]
        bb = gmsh.model.getBoundingBox(2, tag)
        xmin, ymin, zmin, xmax, ymax, zmax = bb
        com = gmsh.model.occ.getCenterOfMass(2, tag)
        cx, cy, cz = com

        # Check adjacent volumes
        # We can simulate this by checking distance to known boundaries
        
        # Inlet (x approx -L_tube)
        if abs(cx - (-L_tube)) < eps:
            inlet_surfs.append(tag)
            continue
            
        # Outlet (x approx L_dom + L_tube)
        if abs(cx - (L_dom + L_tube)) < eps:
            outlet_surfs.append(tag)
            continue
            
        # Heat Source (z approx 0)
        if abs(zmin - 0) < eps and abs(zmax - 0) < eps:
            heat_source_surfs.append(tag)
            continue
            
        # For the rest, we rely on OpenFOAM's splitMeshRegions to find the interface.
        # But we need to define external walls.
        # "walls" usually refers to fluid walls.
        # "solid_walls" refers to solid external.
        
        # Heuristic:
        # If part of fluid volume boundaries AND NOT interface -> Fluid Wall.
        # If part of solid volume boundary AND NOT interface -> Solid Wall.
        
        # This is hard to do purely geometrically without adjacency info.
        # However, OpenFOAM's `autoPatch` or `splitMeshRegions` handles 'defaultFaces' well.
        # We mainly need to tag Inlet, Outlet, and HeatSource because they have special BCs.
        # The rest will default to 'walls' or be split into 'fluid_to_solid'.
        
        # Let's tag 'fluid_walls' as everything else on the outer boundary of fluid domain?
        # Actually, simpler: define patches for specific BCs. The rest becomes defaultFaces (walls).
        pass

    gmsh.model.addPhysicalGroup(2, inlet_surfs, name="inlet")
    gmsh.model.addPhysicalGroup(2, outlet_surfs, name="outlet")
    gmsh.model.addPhysicalGroup(2, heat_source_surfs, name="heat_source")

    # 7. Mesh Generation
    # Refine near walls/interface?
    # Global size set at end
    gmsh.model.mesh.setSize(gmsh.model.getEntities(0), lc_fluid)
    
    gmsh.model.mesh.generate(3)
    
    gmsh.write(output_file)
    gmsh.finalize()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        output_file = sys.argv[1]
    else:
        output_file = "heatsink.msh"
    create_geometry(output_file)
