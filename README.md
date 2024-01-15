# Wearcheck: tribological wear resistance estimator
### Created for "Knowera" by Dmytro Tkachivskyi, PhD
# Table of contents
- [Short description](#Short-description)
- [Screenshots](#Screenshots)
- [Elaborate description](#Elaborate-description)
- [Features](#Features)
    - [Assumptions](#Assumptions)
- [Prerequisites](#Prerequisites)
    - [Getting Python](#Getting-Python)
    - [Getting FreeCAD](#Getting-FreeCAD)
- [Wearcheck package installation](#Wearcheck-package-installation)
- [Adding your own data](#Adding-your-own-data)
    - [Adding new material](#Adding-new-material-materialsjson)
    - [Adding new erodent](#Adding-new-erodent-erodentsjson)
        - [Determining erodent aggressiveness factor](#Determining-erodent-aggressiveness-factor)
- [Story behind Wearcheck](#Story-behind-Wearcheck)
- [Licensing](#Licensing)
- [Support this project](#Support-this-project)

# Short description
**Wearcheck** is a zero‐hassle cross-platform CLI (command line interface) tribological utility for quick abrasive or erosive mechanical wear‐risk estimates on CAD parts. A conjunction of [Archard's equation](https://en.wikipedia.org/wiki/Archard_equation) (modified for erosion and 3-body sliding or abrasion wheel cases), [β-curve distribution function](https://en.wikipedia.org/wiki/Beta_distribution) and [FreeCAD](https://www.freecad.org/) software is used for accurate wear prediction.

# Screenshots

![Erosive wear in terminal](./screenshots/AISI%20316L_erosion_terminal.png)

<p align="center"><i>Image 1. Determining erosion of AISI 316L steel in terminal under 42° impact angle, using theoretical β-distribution curve</i></p>

![Abrasion wear in 3D](./screenshots/Hardox%20400_abrasion.png)

<p align="center"><i>Image 2. Visualising sliding abrasive wear of Hardox 400 steel in 3D</i></p>

![Erosive wear in 3D](./screenshots/TiC-NiMo%20+%20AISI%20316L_erosion.png)

<p align="center"><i>Image 3. Visualising aribitrary heat map during erosive wear of TiC-NiMo + AISI 316L cermet in 3D</i></p>

![Theoretical Erosive wear and heat map in 3D](./screenshots/ICE_parts_erosion.png)

<p align="center"><i>Image 4. Theoretical erosive wear and heat map of ICE's parts with complex geometry. (It would be impossible to check them for the given conditions in real life)</i></p>

# Elaborate description

**Wearcheck** is a single-command “what-if” calculator for mech­anical wear of part of *any* geometrical complexity. It is aimed at engineers who need a *first-order* answer **now**, long before a full DOE or dozens of ASTM/GOST (and other standards) coupons can be scheduled.

Use **Wearcheck** when you *do* care about wear, but *don’t* have the calendar, budget or lab capacity to test every permutation:

* **One command in terminal gives almost instant answer:**  the script converts the STEP to a mesh, shoots rays to find exposed facets, applies an angle-weighted Archard/β-model, and prints:
  * predicted volume loss & risk tier;
  * a plain-text angular wear curve **(abrasion or erosion)**;
  * an optional colour map on the real geometry.
* **Explore the full design space:**  vary load, speed, particle flow, material, impact angle without cutting a single extra coupon.
* **Bridge data gaps:**  when only one ("Time" value for ARWW) or two reference points (30° and 90° for AEW) exist, **Wearcheck** fits a smooth *β-distribution* so you still see the whole curve.
* **No heavyweight CAE or GUI:**  everything is CLI-friendly and runs on commodity hardware; FreeCAD is only called in head-less mode.
* **Fully portable:** the system takes up 2 minutes to set up and works on any machine with Python ≥ 3.9 and FreeCAD ≥ 0.21.

# Features
| Flag&nbsp;/&nbsp;Argument | What it means | Default / Notes |
|---------------------------|---------------|-----------------|
| **Positional STEP / STL / OBJ** | The 3-D part to analyse.  Omit it and the bundled 25 × 50 × 10 mm coupon is used. |
| `--mat, --material` | **Material under test** &nbsp; *(key in `materials.json`)* | mandatory |
| `--wt, --wear-type` | **Mechanism**: `abrasion` (3-body sliding) or `erosion` (particle impact) | abrasion |
| `--impact X Y Z  /  --angle °` | **Attack direction** - either a full 3-D vector (`--impact`) for real parts **or** a simple polar angle (`--angle`) for flat coupons. | required for erosion |
| `--med, --medium` | Abrasive / erodent name &nbsp;*(key in `erodents.json`)* | Silica_Sand |
| `--load` | Normal load **W** (N) – abrasion only | 130 N |
| `--velocity, --vel` | Particle or belt **speed v** (m · s⁻¹) | `v_ref` from erodent DB |
| `--feed-rate, --fr` | Mass flow rate **ṁ** (g · s⁻¹) | mandatory |
| `--time, --t` | Exposure **duration t** (s) | mandatory |
| `--visualize, --vis` | Show a PyVista 3-D scene with risk map, scalar bar, hot/cold toggle, grid & camera shortcuts. | off |
| `--interactive, --int` | Prompt for calibration of custom LOW / MEDIUM / HIGH risk thresholds. | off |
| *Keyboard hot-keys in viewer* | **G** grid, **M** wear ⇄ heat, **I** screenshot, **1–9** preset views. |

## Assumptions

* **Erosive mode**  
  `--angle` ( `-a` ) or `--impact x y z` (`--imp x y z`) describe **only the polar angle θ**
  between the abrasive jet and the horizontal plane.  The auxiliary angles β and δ, and particles flying distance defined in the standard (e.g., GOST 23.201-78) are already absorbed into the fitted Archard/β coefficients; you do **not** enter them separately.

* **Abrasive (sliding) mode**  
  The mechanical wear volume predicted by the 3-body Archard equation is mapped as a *uniform* thinning over the test surface. In real ASTM G65 tests the scar is oval-shaped and centred; the visual map therefore differs, but the **total** wear (mm³) and the wear rate (mm³/N·m) are calibrated to match reference coupon results as closely as possible.  

* **Coated parts (HVOF spraying, cladding, etc.)**  
  When your material is a coating on a substrate (not a solid block like Hardox 400), the predicted wear depth refers to how much of that coating layer is removed. **Always compare** the script’s estimated wear height on the test surface to your known coating thickness — if the wear depth exceeds the layer thickness, the substrate would be exposed in a real test.

# Prerequisites

- Python ≥ 3.9 (to check run "python3 --version" and "pip3 --version");
- FreeCAD ≥ 0.21 (to check run "freecadcmd --version" or "freecadcmd --v).

## Getting Python

- **Windows**  
  Download and run the official installer from https://www.python.org/downloads/windows/ and make sure to check **“Add Python 3.x to PATH”** during installation

- **Linux**  
  ```bash
  sudo apt update && apt install python3 python3-venv python3-pip  # Debian-based

  sudo dnf install python3 python3-virtualenv python3-pip # Fedora/RHEL
  ```
- **MacOS**
  ```bash
  brew install python
  ```

## Getting FreeCAD

- **Windows**  
  Download and run the official installer of FreeCAD desktop package from https://www.freecad.org/downloads.php and make sure to check **“Add FreeCAD to PATH”** during installation
 
- **Linux**
  ```bash
  sudo apt update && sudo apt install freecad # Debian-based
  
  sudo dnf install freecad # Fedora/RHEL
  ```
- **MacOS**
  ```bash
  brew install freecad
  ```

# Wearcheck package installation

If this is familiar for you, clone the repository; otherwise download and unzip the following archive:
```bash
curl -LO https://github.com/dmtkac/knowera-wearcheck/archive/refs/heads/main.zip

# first time using
cd ~/path/to/the/unzipped/archive # switch to where archive was unzipped
python3 -m venv .venv # create virtual environment in the current directory (only first time)
source .venv/bin/activate # activate virtual environment in the current directory
pip install -e . # install dependencies

# subsequent using
cd ~/path/to/the/unzipped/archive # switch to where archive was unzipped
source .venv/bin/activate # activate virtual environment in the current directory

wearcheck --mat "Hardox 400" ... # run with standard 25x50x10 mm step part
wearcheck "/path/to/part.step" --mat "TiC-NiMo + AISI 316L" ... # run with any custom part on the system

wearcheck --help # call help if needed

deactivate # exit virtual environment
```

# Adding your own data 

## Adding new material (`materials.json`)

Fill **7 keys** per material: `hardness_MPa, wear_k, e30, e90, beta_a, beta_b`; *(optionally keep `hardness_HV` or `hardness_HBW` for reference).*

| Step | What you measure / calculate | JSON field(s) you fill |
|------|-----------------------------|------------------------|
| **1.  Hardness** | Any scale is OK (HV, HBW, GPa) &nbsp;→ **convert to MPa**:<br>`MPa = HV × 9.807`  &nbsp;&nbsp;or&nbsp;&nbsp; `MPa = HBW × 9.807` | `"hardness_MPa": …`<br>(optionally, keep the raw number, e.g. `"hardness_HV": …`) |
| **2.  Abrasion coupon** *(ASTM G65)* | Measure coupon volume loss **Q_exp** (mm³).  Note load **W [N]**, sliding distance **L [m]**, abrasive mass **M [kg]**. | |
| **3.  Archard 3-body factor** | `wear_k = (Q_exp · H_MPa) / (W · L · M)` | `"wear_k": …` |
| **4.  Erosion coupons** | Air-jet / slurry tests on flat coupons at **30 °** and **90 °**. Record AEW (mm³ / kg). | `"e30": …`  `"e90": …` |
| **5.  β-curve shape** | Tweak **`beta_a` > 1**, **`beta_b` > 1** until the theoretical curve crosses your e30 & e90 points smoothly. (Use `e30 = e90 = 1` temporarily while tuning.) | `"beta_a": …`  `"beta_b": …` |
| **6.  Insert block** |  e.g., Cr3C2-Ni
```json
"Cr3C2-Ni": {
  "hardness_HV": 1270,
  "hardness_MPa": 12446,
  "wear_k": 41.61e5,
  "e30": 8.956e7,
  "e90": 13.89e7,
  "beta_a": 1.2,
  "beta_b": 1
}
```
**Verify**
```bash
wearcheck "Cr3C2-Ni" --wt erosion   --med "Silica Sand"  --a 42   --vel 40   --fr 2.5   --t 2400
wearcheck "Cr3C2-Ni" --wt abrasion   --med "Silica Sand"  --vel 2.4   --fr 6.17   --t 300
```

## Adding new erodent (`erodents.json`)

Each abrasive / erodent entry needs **only two values**:

| Key        | Meaning                                                         | Example |
|------------|-----------------------------------------------------------------|---------|
| `factor`   | **Relative aggressiveness** (dimensionless).<br>`1.0` is the baseline “Silica Sand”.<br>Wear volume predicted by the model is multiplied by this factor. | Silica Sand → `1.00`<br>White Alumina → `1.73` |
| `v_ref`    | **Reference particle speed** (m · s⁻¹) at which your empirical `e30`/`e90` were measured.<br>The code scales wear with \((v / v_{ref})^{n}\) where *n* defaults to 2 or is given per material. | 40.0 m/s |

### Determining erodent aggressiveness *factor*

1. Run the same coupon test (30 ° or 90 °) with the baseline sand and with your new abrasive;
2. Compute the ratio of volume-loss **per kilogram**: factor = AEW_alumina_30/AEW_silica_30 (or AEW_alumina_90/AEW_silica_90);
3. Round to two decimals and enter as `"factor"`;
4. Add to `erodents.json`:

```json
{
  "Silica Sand"  : { "factor": 1.00, "v_ref": 40.0 },
  "Alumina Wheel": { "factor": 1.73, "v_ref": 60.0 },
  ...
}
```

**How was the Alumina factor determined (1.73)?**

Once Cr3C2–Ni entered the materials library in `materials.json`, it became possible to run wear predictions under identical β-curve calibration for both Silica Sand and an Alumina Wheel at the speed used in AEW [dataset](https://doi.org/10.1016/j.wear.2006.12.027) (in this case 60 m/s) over the 30°–90° range (in this case only AEW_30 and AEW_90 were used). Comparing silica-simulated wear rate to published AEW data, which used alumina wheel, revealed that Cr3C2–Ni lost about 73% more volume per kilogram with alumina than with silica under the same conditions. That directly yields the 1.73 aggressiveness factor for the Alumina Wheel entry.

**Verify**
```bash
wearcheck "Hardox 400" --wt abrasion --med "Alumina Wheel" --vel 3.2 --fr 4.5 --t 600 
wearcheck "TiC-NiMo + AISI 316L" --wt abrasion --med "Alumina Wheel" --vel 2.4 --fr 6.0 --t 300
```

# Story behind Wearcheck

This tool is a product of R&D interest of Knowera enterprise (for whom **Wearcheck** I initially developed it back in 2022) and background from my [PhD research](https://doi.org/10.23658/taltech.41/2021) on sliding and impact wear of HVOF-sprayed hardmetal/steel composites. In 2024 it was decided to share project's version 0.1.0 source code with the scientific/coding community. Most of the **materials database** (Hardox 400, AISI 316L, TiC-NiMo + 316L) comes directly from my own coupon work-ups. An additional entries (material Cr3C2-Ni and erodent Alumina Wheel) were incroprated from my colleagues' research ([I. Hussainova et al., 2007](https://doi.org/10.1016/j.wear.2006.12.027)) to show that the same β-curve/Archard framework adapts to *any* "material + medium" pair once two calibration points are known.

**Wearcheck** is therefore both:

* a **personal lab notebook distilled to code**, and  
* a **generic template** other researchers can extend with their own `materials.json` rows and erodent descriptions.

# Licensing
This work is licensed under the [MIT License](./LICENSE).

# Support this project
If this utility helps you save machine time or sparks ideas for your own experiments, please cite it — or better, contribute your data back to the library or codebase, so it keeps growing, enriching [FOSS](https://en.wikipedia.org/wiki/Free_and_open-source_software) ecosystem.

Also, if you enjoyed this project or found it useful, consider supporting me! Your donations help me to maintain and develop future improvements.

[![Donate with PayPal button](https://www.paypalobjects.com/en_US/i/btn/btn_donate_LG.gif)](https://www.paypal.com/donate?hosted_button_id=WWVH67M22965A)

[![Donate with Crypto](https://img.shields.io/badge/Donate%20with-Crypto-green?style=flat-square&logo=bitcoin)](./qr_crypto.png)