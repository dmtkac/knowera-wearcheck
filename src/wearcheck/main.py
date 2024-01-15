#!/usr/bin/env python3
import scipy, argparse, json, os, shutil, sys, subprocess, tempfile, warnings
import numpy as np, trimesh
from math import gamma, acos, degrees
from pathlib import Path

# FreeCAD check
if shutil.which("freecadcmd") is None:
    sys.exit("FreeCAD not found. Install it first.")

# JSON to Python dict.
def load_json(p):
    with open(p) as f:
        return json.load(f)

# freecadcmd: STEP to STL mesh
def export_step_to_mesh(step):

    if not os.path.isfile(step):
        alt = os.path.join("test_part", step)
        if os.path.isfile(alt):
            step = alt
        else:
            print(f"ERROR: STEP file not found: {step}")
            sys.exit(1)

    # temporary STL
    tf = tempfile.NamedTemporaryFile(suffix=".stl", delete=False)
    tf.close()
    mesh_path = tf.name

    # FreeCAD export
    cmd = f'''
import Part, MeshPart
sh = Part.Shape(); sh.read(r"{step}")
m = MeshPart.meshFromShape(sh, LinearDeflection=0.5, AngularDeflection=1, Relative=False)
m.write(r"{mesh_path}")
'''
    clean_env = os.environ.copy()
    clean_env.pop("VIRTUAL_ENV", None)
    clean_env.pop("PYTHONPATH", None)

    res = subprocess.run(
        ["freecadcmd", "-c", cmd],
        env=clean_env,
        capture_output=True, text=True
    )
    if res.returncode:
        print("Mesh export failed!")
        print("stdout:", res.stdout)
        print("stderr:", res.stderr)
        sys.exit(1)

    return mesh_path

# unit-length vector conversion
def normalize(v):
    n = np.linalg.norm(v)
    return v/n if n>0 else v

# beta-function for theoretical wear distribution
def beta_pdf(x,a=4,b=4):
    B = gamma(a)*gamma(b)/gamma(a+b)
    return x**(a-1)*(1-x)**(b-1)/B

# wear
def main():

    # CLI arguments
    p = argparse.ArgumentParser(
        usage=(
            "wearcheck [FILE] "
            "[--impact X Y Z | --angle ANGLE] "
            "[--material MATERIAL] "
            "[--wear-type {abrasion,erosion}] "
            "[--medium MEDIUM] "
            "[--load LOAD] "
            "[--velocity VELOCITY] "
            "[--feed-rate FEED_RATE] "
            "[--time TIME] "
            "[--visualize] "
            "[--interactive] "
        ),
        description=(
            "Wearcheck: tribological wear resistance estimator (https://github.com/dmtkac/knowera-wearcheck)\n"
            "Created for \"Knowera\" by Dmytro Tkachivskyi, PhD"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    # mutually-excluding --impact and --angle flags  
    impact_grp = p.add_mutually_exclusive_group()
    impact_grp.add_argument("--impact", "--imp", nargs=3, type=float, metavar=("X", "Y", "Z"),
                   help="Impact direction vector [x y z] in part coordinates:\n"
                   " ▸ for real-world parts when the jet can hit several "
                   "faces at different angles (turbines, elbows, impellers, etc.);\n"
                   " ▸ mutually exclusive with \"--angle\"")
    impact_grp.add_argument("--angle", "--a",  type=float,
                   help="Impact angle in degrees from the XY plane "
                   "(0 ° = grazing; 90 ° = normal incidence):\n"
                   " ▸ for ASTM/GOST coupon tests and other single-plane parts;\n"
                   " ▸ mutually exclusive with \"--impact\"")       
    # tribological args 
    p.add_argument("--material", "--mat",  required=True, 
                   help="Material type (specified in \"materials.json\")")
    p.add_argument("--wear-type", "--wt", choices=["abrasion","erosion"], default="abrasion")
    p.add_argument("--medium", "--med",    default="Silica_Sand",
                   help="Abrasive or erodent medium (specified in \"erodents.json\")")
    p.add_argument("--load", "--l",      type=float, default=130.0,
                   help="Normal load (N) for abrasion tests (only for sliding abrasion)"),
    p.add_argument("--velocity", "--vel",  type=float, default=None,
                   help="Relative speed of the abrasive in m/s (if omitted, uses v_ref from \"erodents.json\")")
    p.add_argument("--feed-rate", "--fr", type=float, required=True,
                   help="Abrasive/erosive feed rate in g/s (for abrasive wheel should be set to 0)"),
    p.add_argument("--time", "--t",     type=float, required=True,
                   help="Exposure duration in seconds")
    # visualization args
    p.add_argument("file",  nargs="?", default=None,
                   help="STEP or mesh file; if omitted uses newest file in ./test_part")
    p.add_argument("--visualize", "--vis",  action="store_true",
                   help="Enabling creating of 3D model (optional)")
    p.add_argument("--interactive", "--int", action="store_true",
                   help="Interactive mode for calibrating thresholds (optional)")  
    pos_group = p._action_groups[0]
    opt_group = p._action_groups[1]
    opt_group._group_actions.extend(pos_group._group_actions)
    p._action_groups.remove(pos_group)
    args = p.parse_args()

    # mesh file from 'test_part' dir
    if args.file is None:
        step_dir = Path(__file__).with_name("test_part")
        patterns = ("*.stp", "*.step", "*.stl", "*.obj")
        part_files = [p for pat in patterns for p in step_dir.glob(pat)]
        if not part_files:
            sys.exit("\nERROR: no part files found in test_part/ "
                    "(expected *.step / *.stp / *.stl / *.obj)")
        # priority to the newest file
        newest = max(part_files, key=lambda p: p.stat().st_mtime)
        args.file = str(newest)
        print(f"\nWARNING: using {newest.name} from {step_dir}")

    # flags for special conditions
    if args.wear_type == "erosion" and args.angle is None and args.impact is None:
        p.error("\nERROR: for erosion you must supply either '--impact' or '--angle' flag")
    if args.wear_type == "abrasion":
        print("\nWARNING: '--angle' and '--impact' flags ignored for 'abrasion' mode")
        impact = np.array([0.0, 0.0, -1.0])
        args.angle = 0.0
    else:
        if args.angle is not None:
            θ = np.radians(180-args.angle)
            impact = normalize(np.array([np.cos(θ), 0.0, -np.sin(θ)]))
        else:
            impact = normalize(np.array(args.impact, dtype=float))
    args.impact = impact

    # default thresholds
    low_pct = 5.0
    med_pct = 15.0

    # databases
    mats = load_json("materials.json")
    ers = load_json("erodents.json")
    mat = mats.get(args.material)
    er = ers.get(args.medium)
    if not mat or not er:
        print("Invalid material or medium"); sys.exit(1)

    # beta_a and beta_b validation
    a = mat.get("beta_a")
    b = mat.get("beta_b")
    if a < 1.0 or b < 1.0:
        msg = (f"\nERROR: beta_a ({a}) and beta_b ({b}) must both be >= 1 in 'materials.json' "
            f"for a physically valid wear-angle distribution.\n"
            f"Please refit the curve with β-parameters >= 1.")
        sys.exit(msg)

    # step to mesh
    ext = os.path.splitext(args.file)[1].lower()
    if ext in (".step",".stp"):
        mesh_file = export_step_to_mesh(args.file)
    else:
        mesh_file = args.file

    # part volume
    mesh = trimesh.load(mesh_file, force="mesh")
    min_face_area = 1e-6
    if (mesh.area_faces < min_face_area).any():
        keep = mesh.area_faces >= min_face_area
        mesh.update_faces(keep)
        mesh.remove_unreferenced_vertices()
    if mesh.is_empty:
        print("Failed to load mesh"); sys.exit(1)
    if mesh.is_watertight:
        vol = mesh.volume
    else:
        vol = mesh.convex_hull.volume
        print(f"Warning: using convex hull volume {vol:.2f} mm³")

    # parameters in terminal
    print()
    print(f"Part name: {Path(args.file).name}")
    print(f"Material: {args.material}")
    print(f"Wear type: {args.wear_type}")
    if args.wear_type == "abrasion":
        print(f"Impact vector: N/A (sliding)")    
    else:
        print(f"Impact vector: {normalize(np.array(args.impact))}")
    if args.wear_type == 'abrasion':
        print(f"Duration: {args.time} s")
    print(f"Total part volume: {vol:.2f} mm³")

    # ray per 1 mesh triangle
    impact = normalize(np.array(args.impact))
    R = mesh.extents.max()*2.0
    cents = mesh.triangles_center
    origins = cents - impact*R
    dirs = np.tile(impact, (len(cents),1))
    first = mesh.ray.intersects_first(origins, dirs)

    # exposed mesh triangles
    exposed = []
    for i, tri in enumerate(mesh.triangles):
        if first[i] == i:
            area = mesh.area_faces[i]
            n = normalize(mesh.face_normals[i])
            angle = np.degrees(np.arccos(np.clip(n.dot(impact), -1,1)))
            exposed.append((area, angle, tri))

    # wear loss
    if args.velocity is None:
        if "v_ref" in er:
            args.velocity = er["v_ref"]
            print(f"\nUsing v_ref from erodents.json: {args.velocity:.1f} m/s")
        else:
            args.velocity = 40.0
            print(f"\nNo velocity supplied or in erodents.json; defaulting to {args.velocity:.1f} m/s")
    else:
        print(f"Particles velocity: {args.velocity:.1f} m/s")

    # edge case for abrasive wheel
    is_wheel = False 
    if args.feed_rate == 0:
        total_mass = 1.0
        is_wheel = True
        print("Abrasive wheel selected (feed-rate 0 g/s); "
          "mass term normalised to 1.")
    else:
        total_mass = args.feed_rate / 1000.0 * args.time   
        print(f"Abrasive mass used: {total_mass:.3f} kg "
            f"(feed-rate {args.feed_rate} g/s over {args.time}s)")


    # empirical wear metrics
    AEW = None # erosion metric (mm³ / kg)
    ARWW = None # abrasion metric (mm³ / N·m)
    
    # Vickers hardness as pressure (Pa)
    H = mat["hardness_MPa"] * 1e6

    # total abrasive mass available (kg)
    M = total_mass

    if args.wear_type == "abrasion":

        # total sliding distance (m)
        L = args.velocity * args.time
        
        # normal load (N)
        W = args.load

        # Archard's coefficient (obtained empirically from 3-body tests)
        # used for calibrating the model
        K3b = mat["wear_k"]

        # Archard's abrasion volume (adjusted for 3-body and abrasion wheel tests)
        Q = K3b * W * L / H
        if not is_wheel:
            Q *= total_mass

        # scale factor of abrasive media (in erodents.json)
        E_scale = er.get("factor", 1)

        # final abrasion volume (adjusted for different abrasive media)
        score = Q * E_scale
        ARWW = score / (L * W)
        print(f"\nCalculated ARWW value: {ARWW:.2e} mm³/N·m")

        # for 3D visual (per-triangle wear, mm³):
        cell_wear = np.zeros(len(mesh.faces), dtype=float)
        
    else:

        # experimental constants and global scale factors
        v_ref = er.get("v_ref", args.velocity)
        v_expo = er.get("velocity_exponent", 2)
        v_scale = (args.velocity / v_ref) ** v_expo
        media_factor= er.get("factor", 1)
        
        # faces actually hit by the abrasive
        exp_idx = np.array([i for i,hit in enumerate(first) if hit == i])
        
        # geometry of exposed faces
        n_exp = mesh.face_normals[exp_idx]
        areas_exp = mesh.area_faces[exp_idx]
        
        # incidence angle and material response
        cos_theta = n_exp @ impact
        surf_inc = np.abs(90.0 - np.degrees(np.arccos(np.clip(cos_theta, -1, 1))))
        x = surf_inc / 90.0 
        a, b = mat["beta_a"], mat["beta_b"]
        e30, e90 = mat["e30"], mat["e90"]
        w = beta_pdf(x, mat["beta_a"], mat["beta_b"])
        w_max = w.max()
        E = w * e30 + (1.0 - w) * e90

        # contribution of each triangle to overall erosion
        contrib = areas_exp * E * cos_theta**2 / H 
        AEW = contrib.sum()

        # actual volume loss per face (mm³)
        tri_wear = contrib * M * v_scale * media_factor

        # per-face wear in array
        cell_wear = np.zeros(len(mesh.faces))
        cell_wear[exp_idx] = tri_wear

        # total predicted volume loss for the whole part (adjusted for erodent velocity and impact angles)
        AEW = max(AEW, 0.0)
        if args.angle:
            print(f"\nImpact angle: {args.angle}\nCalculated AEW value: {AEW:.2f} mm³/kg")
        else:
            print(f"\nCalculated AEW value: {AEW:.2f} mm³/kg")
        score = AEW * M * v_scale * media_factor

    # percentage computation and classification
    pct = (score / vol * 100) if vol > 0 else 0

    # custom thresholds (if needed)
    if args.interactive:
        val = input(f"LOW threshold (%) [{low_pct}]: ")
        low_pct = float(val) if val.strip() else low_pct
        val = input(f"MEDIUM threshold (%) [{med_pct}]: ")
        med_pct = float(val) if val.strip() else med_pct
        print(f"Using thresholds LOW={low_pct}%, MEDIUM={med_pct}%\n")

    # colour by risk
    BOLD = "\033[1m"
    RED = "\033[31m"
    RESET = "\033[0m"

    if pct < low_pct:
        lvl = 'LOW'
    elif pct < med_pct:
        lvl = 'MEDIUM'
    else:
        lvl = 'HIGH'

    if lvl == "LOW":
        color = "\033[32m"   # green
    elif lvl == "MEDIUM":
        color = "\033[33m"   # yellow
    else:
        color = "\033[31m"   # red

    print(f"{BOLD}{color}Predicted wear: {score:.2e} mm³ ({pct:.2f}%)  |  Risk: {lvl}{RESET}")

    # ASCII distribution curve (in terminal)
    if args.wear_type == "erosion" and args.angle is not None:
        angles = np.linspace(0, 90, 200)
        x_plot = angles / 90.0
        w_plot = beta_pdf(x_plot, mat["beta_a"], mat["beta_b"])
        w_plot /= w_plot.max()
        E_plot = w_plot * e30 + (1.0 - w_plot) * e90
        width, height = 60, 10
        angles_plot = np.linspace(0, 90, width)
        global_AEW = []
        R = mesh.extents.max()*2.0
        centers = mesh.triangles_center
        M_factor = (args.feed_rate/1000*args.time) * ( (args.velocity/er.get("v_ref",args.velocity))**er.get("velocity_exponent",2) ) * er.get("factor",1)
        for θ in angles_plot:
            th = np.radians(180-θ)
            imp = normalize(np.array([np.cos(th), 0.0, -np.sin(th)]))
            origins = centers - imp*R
            first = mesh.ray.intersects_first(origins, np.tile(imp, (len(centers),1)))
            hits = np.flatnonzero(first == np.arange(len(first)))
            n_exp = mesh.face_normals[hits]
            areas_exp = mesh.area_faces[hits]
            cos_t = n_exp.dot(imp)
            surf_inc = np.abs(90 - np.degrees(np.arccos(np.clip(cos_t, -1,1))))
            x = surf_inc/90.0
            w = beta_pdf(x, a, b)
            E = w*e30 + (1-w)*e90
            G = (areas_exp * E * cos_t**2 / H).sum()
            global_AEW.append(max(G,0.0))
        ymin, ymax = np.min(global_AEW), np.max(global_AEW)
        rows = [[' ']*width for _ in range(height)]
        pos30 = int(30/90*(width-1))
        pos90 = width - 1
        pos_cur = None 
        if args.angle is not None:
            angle_cur = args.angle
        else:
            angle_cur = np.degrees(np.arcsin(abs(impact[2])))
        pos_cur = int(angle_cur/90*(width-1))
        prev_r = None
        for col, val in enumerate(global_AEW):
            r = height - 1 - int(round((val - ymin) / (ymax - ymin + 1e-12) * (height - 1)))
            if prev_r is not None:
                step = 1 if r > prev_r else -1
                for rr in range(prev_r, r, step):
                    if rows[rr][col] == ' ':
                        rows[rr][col] = '•'
            if col == pos_cur:
                rows[r][col] = '◆'
            elif col in (pos30, pos90):
                rows[r][col] = '▼'
            else:
                rows[r][col] = '•'
            prev_r = r
        label_w = max(6, len(f"{ymax:.2f}"))
        tick_rows = {}
        for i in range(10):
            r = int(i*(height-1)/9)
            val = ymax - (ymax - ymin)*i/9
            tick_rows[r] = f"{val:>{label_w}.2f} ┤ "
        blank_prefix = " " * label_w + " │ "
        print(f"\nSpecific AEW distribution curve for {args.material} (only for flat surfaces; standard: GOST 23.201-78):\n")
        print("AEW, mm³/kg")
        print(" " * label_w + " ^")
        for r in range(height):
            prefix = tick_rows.get(r, blank_prefix)
            row_str = "".join(rows[r])
            shifted = row_str[2:] + " "
            print(prefix + shifted)
        print(" " * label_w + " └" + "─"*width + " > θ,°")
        xticks = [" "] * width
        for angle in range(0, 91, 15):
            pos = int(angle/90*(width-1))
            xticks[pos] = "|"
        print(" " * (label_w + 1) + "".join(xticks))
        xlabel = [" "] * width
        for angle in range(0, 91, 15):
            lbl = f"{angle:2d}"
            pos = int(angle/90*(width-1)) - len(lbl)//2
            for i,ch in enumerate(lbl):
                if 0 <= pos+i < width:
                    xlabel[pos+i] = ch
        print(" " * (label_w + 1) + "".join(xlabel))
        if args.angle:
            print(f"\n   •  theoretical curve  |  ◆  predicted AEW ({args.angle}°)  |  ▼  empirical AEW (30°, 90°)\n")
        else:
            v = impact
            abs_v = np.abs(v)
            max_c = abs_v.max()
            tol = 1e-6
            planes = [f"{np.degrees(np.arcsin(c)):.1f}°"
                    for i, c in enumerate(abs_v) if c >= max_c - tol]
            plane_str = ", ".join(planes)
            print(f"\n   •  theoretical curve  |  ◆  predicted AEW ({plane_str})  |  ▼  empirical AEW (30°, 90°)\n") 
    elif args.wear_type == "abrasion": 
        width = 60
        height = 15
        t_final = args.time
        ref_time = args.time + 60.0 
        t_max = max(t_final, ref_time)
        slope = ARWW / t_final
        arww_max = slope * t_max
        times = np.linspace(0, t_max, width)
        arww = slope * times 
        epsilon = 1e-12
        step = max(arww_max, epsilon) / (height - 1)
        rows = [[' ']*width for _ in range(height)]
        for col, val in enumerate(arww):
            r = height - 1 - int(val / step + 1e-12)
            rows[r][col] = '•'
        col_star = int(t_final / t_max * (width - 1))
        r_star = height - 1 - int(ARWW / step + 1e-12)
        rows[r_star][col_star] = '▼'
        label_w = max(6, len(f"{arww_max:.2e}"))
        tick = [f"{(i*step):>{label_w}.2e} ┤ " for i in range(height)]
        tick.reverse()
        if args.feed_rate ==0:
            print(f"\nExtrapolated ARWW-Time dependency for {args.material} (standard: Abrasive wheel used):\n")
        else:
            print(f"\nExtrapolated ARWW-Time dependency for {args.material} (standard: ASTM G65):\n")
        print("ARWW, mm³/N·m")
        print(" " * label_w + " ^")
        for r in range(height):
                row_str = "".join(rows[r])
                shifted = row_str[2:] + " "
                print(tick[r] + shifted if r < len(tick) else blank + shifted)
        print(" " * label_w + " └" + "─"*width + " > Time, s")
        def put_label(fraction: float, text: str, buffer: list[str]):
            col = int(fraction*(width-1)) - len(text)//2
            col = max(0, min(col, width - len(text)))
            for i, ch in enumerate(text):
                buffer[col + i] = ch
        xticks = [" "] * width
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            put_label(frac, "|", xticks)
        print(" " * (label_w + 1) + "".join(xticks))
        labels = [" "] * width
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            put_label(frac, f"{int(frac*t_max):d}", labels)
        print(" " * (label_w + 1) + "".join(labels))
        print(f"\n   •  predicted ARWW  |  ▼  empirical ARWW at {t_final}s\n")

# CAD-style 3D model (separate window)
    if args.visualize:
        try:
            import pyvista as pv
        except ModuleNotFoundError:
            print("PyVista not installed. Run 'pip install pyvista' in your venv, or omit --visualize to skip viewing.")
            return
        try:
            # PolyData
            faces = np.hstack((np.full((mesh.faces.shape[0],1), 3, dtype=np.int64), mesh.faces.astype(np.int64)))
            pv_mesh = pv.PolyData(mesh.vertices, faces)
            pv_mesh.cell_data['wear'] = cell_wear

            # plotter window
            if args.wear_type == "abrasion":
                impact_str = "N/A (sliding)"
            else:
                impact_str = " ".join(f"{v:.3f}" for v in args.impact)
            title_str = (
                f"Part name: {Path(args.file).name} | "
                f"Part material: {args.material} | "
                f"Medium: {args.medium} | "
                f"Wear type: {args.wear_type} | "
                f"Impact vector: {impact_str}"
            )
            plotter = pv.Plotter(window_size=(1366,768), title=title_str)
            plotter.set_background("whitesmoke", top="steelblue")
            plotter.add_axes()  
            # legend
            if args.wear_type == "abrasion":
                plane_str = "XY at 0.0°"
            elif args.angle is not None:
                plane_str = f"XY at {args.angle:.1f}°"
            else:
                v = impact
                abs_v = np.abs(v)
                max_c = abs_v.max()
                tol = 1e-6
                planes = [("YZ","ZX","XY")[i] + f" at {np.degrees(np.arcsin(c)):.1f}°"
                        for i, c in enumerate(abs_v) if c >= max_c - tol]
                plane_str = ", ".join(planes)
            entries = [
                ("Predicted wear", f"{score:.2e} mm³ ({pct:.2f}%)"),
            ]
            if args.wear_type == "abrasion":
                entries.append(("Calculated ARWW", f"{ARWW:.2e} mm³/N·m")),
                if args.feed_rate == 0: entries.append(("Standard",  "Abrasive wheel used")) 
                else: entries.append(("Standard",  "ASTM G65"))
            else:
                entries.append(("Calculated AEW",  f"{AEW:.2f} mm³/kg")),
                entries.append(("Standard",  "GOST 23.201-78"))
            entries.append(("Part name", Path(args.file).name))
            entries.extend([
                ("Tested material",  args.material),
                ("Part volume",  f"{vol:.2f} mm³"),
                ("Abrasive mater.",  args.medium),
            ])
            if args.wear_type == "abrasion":
                entries.append(("Angle with hor.",  f"N/A (sliding)"))
            else:
                entries.append(("Angle with hor.",  plane_str))
            entries.extend([
                ("Abrasive mass", f"{total_mass:.3f} kg"),
                ("Particle vel.", f"{args.velocity:.1f} m/s"),
                ("Wear type",  args.wear_type),
                ("Risk",  lvl),
            ])
            max_lbl = max(len(lbl) for lbl,_ in entries)
            lines = [
                f"{lbl.ljust(max_lbl)} : {val}"
                for lbl, val in entries
            ]
            info = "\n".join(lines)
            text_actor = plotter.add_text(
                info,
                position=(0.01, 0.985),
                font_size=11,
                color='white',
                shadow=True,
                font='courier',
                viewport=True
            )
            tp = text_actor.GetTextProperty()
            tp.SetJustificationToLeft() 
            tp.SetVerticalJustificationToTop() 

            # colormap distribution
            if args.wear_type == "abrasion":
                nz = np.abs(mesh.face_normals[:, 2])
                xy_faces = np.where(nz > 0.95)[0]
                total_xy_area = mesh.area_faces[xy_faces].sum()
                for fi in xy_faces:
                    share = mesh.area_faces[fi] / total_xy_area
                    cell_wear[fi] = score * share
                abrasion_thickness = cell_wear / mesh.area_faces
                thickness = abrasion_thickness
                vertex_thickness = np.zeros(mesh.vertices.shape[0])
                counts = np.zeros_like(vertex_thickness)
                for fi in xy_faces:
                    for vid in mesh.faces[fi]:
                        vertex_thickness[vid] += abrasion_thickness[fi]
                        counts[vid] += 1
                mask = counts > 0
                vertex_thickness[mask] /= counts[mask]
            else: # (erosion)
                if args.angle is not None:
                    nz = np.abs(mesh.face_normals[:, 2])
                    xy_mask = nz > 0.95
                    cell_wear[~xy_mask] = 0.0
                thickness = cell_wear / mesh.area_faces
                warp_factor = -0.01
            n_pts = mesh.vertices.shape[0]
            vertex_thickness = np.zeros(n_pts)
            counts = np.zeros(n_pts)
            exposed_idxs = [i for i,hit in enumerate(first) if hit == i]
            for fi in exposed_idxs:
                for vid in mesh.faces[fi]:
                    vertex_thickness[vid] += thickness[fi]
                    counts[vid] += 1
            mask = counts > 0
            vertex_thickness[mask] /= counts[mask]
            impact_dir = impact
            if args.angle is not None:
                tol = 1e-3
                primary = [i for i,n in enumerate(mesh.face_normals) if abs(abs(n[2]) - 1.0) < tol]
            else:
                dp_all = mesh.face_normals.dot(impact_dir)
                best_val = dp_all.min()
                tol = 1e-3
                primary = [i for i,v in enumerate(dp_all) if abs(v - best_val) < tol]
            if len(primary) == 0:
                primary = exposed_idxs
            areas = mesh.area_faces[primary]
            centers = mesh.triangles_center[primary]
            center_pt = (centers * areas[:,None]).sum(axis=0) / areas.sum()
            pts = mesh.vertices
            dP = pts - center_pt
            proj = dP - np.outer(dP.dot(impact_dir), impact_dir)
            radial = np.linalg.norm(proj, axis=1)
            
            # arrows
            spread = mesh.extents.max() * 0.15
            hover = mesh.extents.max() * 0.6
            length = mesh.extents.max() * 0.2
            if args.wear_type == "abrasion":
                z_surface = mesh.bounds[1][2]
                y_min, y_max = mesh.bounds[0][1], mesh.bounds[1][1]
                y_positions = np.linspace(y_min + 0.25*(y_max-y_min), y_max - 0.25*(y_max-y_min), 3)
                x_tail = -mesh.bounds[0][0] + length
                rotated_arrows = np.array([impact[2], impact[1], impact[0]])
                for y in y_positions:
                    tail = np.array([x_tail, y, z_surface])
                    arrow = pv.Arrow(start=tail, direction=rotated_arrows, scale=length)
                    plotter.add_mesh(arrow, color="purple")
            else: # erosion
                up = np.array([0, 0, 1])
                if abs(impact.dot(up)) > 0.9:
                    up = np.array([0, 1, 0])
                u = normalize(np.cross(impact, up))
                v = normalize(np.cross(impact, u))
                for dx, dy in [( 1,  0), (-1,  0), ( 0, 0)]:
                    tip = center_pt + u*(dx*spread) + v*(dy*spread)
                    tail = tip - impact_dir*hover
                    arrow = pv.Arrow(start=tail, direction=impact_dir, scale=length)
                    plotter.add_mesh(arrow, color='purple')

            # scale calibration
            # wear depth
            if args.wear_type == "abrasion":
                beam_thickness = vertex_thickness * 2
            else:
                p = 1.5 if args.angle is not None else 0.01
                beam_thickness = vertex_thickness * (radial / radial.max())**p
            pv_mesh.point_data['beam_thickness'] = beam_thickness
            # arbitrary surface heat         
            surface_heat = beam_thickness.copy()
            pv_mesh.point_data['surface_heat'] = surface_heat

            # colormap
            pv_mesh = pv_mesh.subdivide(2, 'linear')
            if args.wear_type == "abrasion":
                warped = pv_mesh
            else:
                warped = pv_mesh.warp_by_scalar('beam_thickness', factor=warp_factor)
            # colormap and scale switch
            cmap_cycle = ['cividis', 'turbo']
            current_idx = [0] 
            cold_label = None
            warm_label = None
            hot_label  = None
            def helper(name):
                return pv.LookupTable(cmap=name)
            def toggle_cmap():
                nonlocal cold_label, warm_label, hot_label
                current_idx[0] = 1 - current_idx[0]
                if cold_label:
                    plotter.remove_actor(cold_label)
                    cold_label = None
                if warm_label:
                    plotter.remove_actor(warm_label)
                    warm_label = None
                if hot_label:
                    plotter.remove_actor(hot_label)
                    hot_label = None
                sb = plotter.scalar_bar
                if current_idx[0] == 0:
                    lut = helper('cividis')
                    actor.mapper.SetArrayName('beam_thickness')
                    sb.SetTitle('Wear depth (mm)\n')
                    sb.SetNumberOfLabels(4)
                    tp = sb.GetLabelTextProperty()
                    tp.SetOpacity(1.0)
                    tp.SetFontSize(20)
                    sb.SetDrawAnnotations(False)
                else:
                    lut = helper('turbo')
                    actor.mapper.SetArrayName('surface_heat')
                    sb.SetTitle('Surface t. (au)\n')
                    sb.SetNumberOfLabels(1)
                    tp = sb.GetLabelTextProperty()
                    tp.SetOpacity(0.0)
                    sb.SetDrawAnnotations(False)
                    cold_label = plotter.add_text(
                        "Ambient",
                        position=(0.901, 0.613),
                        font_size=11,
                        color='white',
                        font='courier',
                        shadow=True,
                        viewport=True
                    )
                    warm_label = plotter.add_text(
                        "Warm",
                        position=(0.901, 0.613 + 0.1415),
                        font_size=11,
                        color='white',
                        font='courier',
                        shadow=True,
                        viewport=True
                    )
                    hot_label = plotter.add_text(
                        "Hot",
                        position=(0.901, 0.613 + 0.283), 
                        font_size=11,
                        color='white',
                        font='courier',
                        shadow=True,
                        viewport=True
                    )
                actor.mapper.lookup_table = lut
                sb.SetLookupTable(lut)
                plotter.render()
            plotter.add_key_event('m', toggle_cmap)
            actor = plotter.add_mesh(
                warped,
                scalars='beam_thickness',
                cmap=cmap_cycle[current_idx[0]],
                show_edges=False,
                lighting=True,
                smooth_shading=True,
                interpolate_before_map=False,
                scalar_bar_args={
                    'title': 'Wear depth (mm)\n',
                    'n_labels': 4,
                    'fmt': '%.2e',
                    'vertical': True,
                    'position_x': 0.88,
                    'position_y': 0.62,
                    'width': 0.05,
                    'height': 0.35,
                    'color': 'white',
                }
            )
            # extras
            # custom font
            from vtk import vtkScalarBarActor
            sb_actor = None
            for a in plotter.renderer.GetActors2D():
                if isinstance(a, vtkScalarBarActor):
                    sb_actor = a
                    break
            assert sb_actor, "Could not find the scalar bar actor!"
            title_tp = sb_actor.GetTitleTextProperty()
            label_tp = sb_actor.GetLabelTextProperty()
            title_tp.SetFontFamilyToCourier()
            label_tp.SetFontFamilyToCourier()
            title_tp.SetFontSize(26)
            label_tp.SetFontSize(20)
            title_tp.ShadowOn()
            label_tp.ShadowOn()
            plotter.render()
            # wireframe with part's dimensions
            xmin, xmax, ymin, ymax, zmin, zmax = warped.bounds
            dx = xmax - xmin 
            dy = ymax - ymin
            dz = zmax - zmin
            bbox = pv.Box(bounds=(xmin, xmax, ymin, ymax, zmin, zmax))
            plotter.add_mesh(bbox, color='black', style='wireframe', line_width=2)
            labels = {
                f"X = {dx:.1f} mm": np.array([(xmin+xmax)/-2, ymin, zmin]),
                f"Y = {dy:.1f} mm": np.array([xmin, (ymin+ymax)/-2, zmin]),
                f"Z = {dz:.1f} mm": np.array([xmin, ymin, (zmin+zmax)/-2]),
            }
            pts = -np.vstack(list(labels.values()))
            txts = list(labels.keys())
            plotter.add_point_labels(
                pts,
                txts,
                font_size=16,
                text_color='black',
                point_size=0,
                shape='rounded_rect',
                shape_color='white',
                shape_opacity=0.5
            )
            # grid
            plotter.camera.Zoom(1.1)
            xmin, ymin, zmin = mesh.bounds[0]
            xmax, ymax, zmax = mesh.bounds[1]
            grid_z = zmin - 0.1
            margin = 0.3
            n_target = 12 
            side_xy = max(xmax - xmin, ymax - ymin)
            grid_side = side_xy * (1 + 2 * margin)
            plane = pv.Plane(
                center = ((xmin + xmax) / 2, (ymin + ymax) / 2, grid_z),
                direction = (0, 0, 1),
                i_size = grid_side,
                j_size = grid_side,
                i_resolution = n_target,
                j_resolution = n_target
            )
            grid_actor = plotter.add_mesh(plane, style = 'wireframe', color = 'gray', line_width = 0.5)
            # grid on/off switch
            def toggle_grid():
                vis = grid_actor.GetVisibility()
                grid_actor.SetVisibility(not vis)
                plotter.add_axes()
                plotter.render()
            plotter.add_key_event('g', toggle_grid)
            # screenshot
            import datetime
            def take_screenshot(_=None):
                ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                fname = f"{args.material}_{args.wear_type}_{ts}.png"
                plotter.screenshot(fname)
                toast = plotter.add_text(
                    "Image saved",
                    position=(0.5, 0.87),
                    font_size=14,
                    color='darkgreen',
                    shadow=False,
                    font='courier',
                    viewport=True
                )
                vtk_iren = plotter.ren_win.GetInteractor()
                def clear_toast(caller, event):
                    plotter.remove_actor(toast, reset_camera=False)
                    plotter.add_axes()
                    plotter.render() 
                    vtk_iren.RemoveObserver(tag)
                tag = vtk_iren.AddObserver("TimerEvent", clear_toast)
                vtk_iren.CreateOneShotTimer(1000)
                print(f"Saved image: {fname}")
            plotter.add_key_event('i', take_screenshot) 
            # key-buttons
            def with_grid_hidden(func):
                def wrapped():
                    was_visible = grid_actor.GetVisibility()
                    grid_actor.SetVisibility(False)
                    func()
                    grid_actor.SetVisibility(was_visible)
                    plotter.render()
                return wrapped
            @with_grid_hidden 
            def view_back():
                plotter.view_vector((-1, 0, 0), viewup=(0, 0, 1))
                plotter.camera.Zoom(1.1)
            @with_grid_hidden 
            def view_top():
                plotter.view_vector((0, 0, 1), viewup=(0, -1, 0))
                plotter.camera.Zoom(1.1)
            @with_grid_hidden
            def view_right():
                plotter.view_vector((0, 1, 0), viewup=(0, 0, 1))
                plotter.camera.Zoom(1.1)
            @with_grid_hidden
            def view_left():
                plotter.view_vector((0, -1, 0), viewup=(0, 0, 1))
                plotter.camera.Zoom(1.1)
            @with_grid_hidden
            def view_front():
                plotter.view_vector((1, 0, 0), viewup=(0, 0, 1))
                plotter.camera.Zoom(1.1)
            @with_grid_hidden
            def view_bottom():
                plotter.view_vector((0, 0, -1), viewup=(0, -1, 0))
                plotter.camera.Zoom(1.1)
            @with_grid_hidden
            def view_isom():
                plotter.view_isometric()
                plotter.camera.Zoom(1.1)
                plotter.camera.SetParallelProjection(False)
            keymap = {
                "9": view_back,  "8": view_top,    "6": view_right,
                "5": view_front, "4": view_left,   "2": view_bottom,
                "1": view_isom,
            }
            for k, func in keymap.items():
                plotter.add_key_event(k,       func)
                plotter.add_key_event(f"KP_{k}", func)
            entries = [
                ("G", "toggle grid"),
                ("M", "toggle cmap"),
                ("I", "save image"),
                ("9", "back view"),
                ("8", "top view"),
                ("6", "right view"),
                ("5", "front view"),
                ("4", "left view"),
                ("2", "bottom view"),
                ("1", "isom. view"),
            ]
            max_lbl = max(len(lbl) for lbl,_ in entries)
            lines = [
                f"{lbl.ljust(max_lbl)} : {val}"
                for lbl, val in entries
            ]
            info = "\n".join(lines)
            text_actor = plotter.add_text(
                info,
                position=(0.855, 0.01),
                font_size=10,
                color='black',
                viewport=True ,
                shadow=True,
                font='courier'
            )
            tp = text_actor.GetTextProperty()
            tp.SetJustificationToLeft()
            tp.SetVerticalJustificationToBottom()
            # disclaimer (for abrasion only)
            if args.wear_type == "erosion" and args.angle is None:
                disc_msg = (
                    "Experimental mode for complex parts:\n"
                    "colormaps may be not 100% accurate"
                )
                plotter.add_text(
                    disc_msg,
                    position=(0.45, 0.94),
                    font_size=9,
                    color="darkred",
                    shadow=False,
                    viewport=True,
                    font="courier"
                )
            elif args.wear_type == "abrasion":
                disc_msg = (
                    "Homogeneous wear \n"
                    "distrib. assumed"
                )
                plotter.add_text(
                    disc_msg,
                    position=(0.5, 0.94),
                    font_size=9,
                    color="darkred",
                    shadow=False,
                    viewport=True,
                    font="courier"
                )
            plotter.show()

        except Exception as e:
            print("Visualization failed:", e)

    if ext in (".step",".stp"):
        os.remove(mesh_file)

if __name__ == "__main__":
    main()