import numpy as np

vanjski_promjer = 160   # mm
promjer_grla_mlaznice = 2 * 16.7   # mm

tlak_komore = 30 #bar
tlak_fluida = 70 #bar
cvrstoca_materijala = 230 #MPa
faktor_sigurnosti = 1.5

debljina_vanjskog_zida_ulaz = round((faktor_sigurnosti * tlak_fluida* 1e5 * vanjski_promjer / 1e3) / (2*cvrstoca_materijala * 1e6),4) #m

d_hidraulicke_cijevi = 25    # mm
maseni_protok_IPA = 0.29    # kg/s
gustoca_IPA = 785   # kg/m3
n = 12
k = 1

debljina_Stjenke_Rebra_NaGrlu_Mlaznice = 2 #mm


povrsina_cijevi = np.pi * (d_hidraulicke_cijevi/1000)**2 /4
brzina_protoka_IPA = maseni_protok_IPA/(gustoca_IPA*povrsina_cijevi)
print(f"brzina_protoka_IPA-cijev: {brzina_protoka_IPA:.2f} m/s\n")
Auk = ((d_hidraulicke_cijevi/1000)**2 /4)*np.pi                                #m^2

          
min_cirkularna_debljina_kanala = np.sqrt(Auk*4/np.pi + (promjer_grla_mlaznice/1000)**2) * 1000 - promjer_grla_mlaznice #debljina kanala pri kojoj nema ni jedan zarez aluminijski

debljina_unutarnjeg_zida_grlo = max(round(((faktor_sigurnosti*1.2) * (tlak_fluida - tlak_komore)* 1e5 * promjer_grla_mlaznice/1e3) / (2*cvrstoca_materijala * 1e6),4),1.2 * 1e-3) #m
print(f"debljina_unutarnjeg_zida_grlo: {debljina_unutarnjeg_zida_grlo*1e3} mm")
unutarnja_stjenka_rashladnih_kanala = promjer_grla_mlaznice + 2 * debljina_unutarnjeg_zida_grlo * 1e3

fi_alu_rad = 2 * (debljina_Stjenke_Rebra_NaGrlu_Mlaznice) / (unutarnja_stjenka_rashladnih_kanala)
fi_alu_stupnjevi = np.degrees(fi_alu_rad)
print(f"fi_alu_rad: {fi_alu_rad:.6f} rad")
print(f"fi_alu_stupnjevi: {fi_alu_stupnjevi:.6f}°\n")


d_vani_grlo = np.sqrt((Auk * (1+k) * 4/np.pi + ((promjer_grla_mlaznice/1000)+2*debljina_unutarnjeg_zida_grlo)**2 * (1 - n *fi_alu_rad/(np.pi)))/(1 - n *fi_alu_rad/(np.pi))) # m
print(f"d_vani_grlo: {d_vani_grlo*1000:.2f} mm")
cirkularna_debljina_kanala = (d_vani_grlo - promjer_grla_mlaznice/1000)/2
print(f"cirkularna_debljina_kanala: {cirkularna_debljina_kanala*1000:.2f} mm")

debljina_vanjskog_zida_grlo = round(((faktor_sigurnosti*1.1) * tlak_fluida* 1e5 * d_vani_grlo) / (2*cvrstoca_materijala * 1e6),4) #m
print(f"debljina_vanjskog_zida_grlo: {debljina_vanjskog_zida_grlo*1e3} mm\n")


fi_IPA_ulaz_rad = Auk * (8/n) * (1/(d_vani_grlo**2 - (unutarnja_stjenka_rashladnih_kanala/1000)**2))
fi_IPA_ulaz_deg = np.degrees(fi_IPA_ulaz_rad)
duzina_IPA_ulaz_kanala = fi_IPA_ulaz_rad * promjer_grla_mlaznice / 2

print(f"fi_IPA_ulaz_deg: {fi_IPA_ulaz_deg:.6f}°")
print(f"duzina_IPA_ulaz_kanala: {duzina_IPA_ulaz_kanala:.2f} mm")

fi_IPA_izlaz_rad = k * Auk * (8/n) * (1/(d_vani_grlo**2 - (unutarnja_stjenka_rashladnih_kanala/1000)**2))
fi_IPA_izlaz_deg = np.degrees(fi_IPA_izlaz_rad)
duzina_IPA_izlaz_kanala = fi_IPA_izlaz_rad * promjer_grla_mlaznice / 2

print(f"fi_IPA_izlaz_deg: {fi_IPA_izlaz_deg:.6f}°")
print(f"duzina_IPA_izlaz_kanala: {duzina_IPA_izlaz_kanala:.2f} mm\n")

test360 = n *(2 * fi_alu_stupnjevi + fi_IPA_ulaz_deg + fi_IPA_izlaz_deg)
print(f"test-grlo-360: {test360:.6f}°\n")


# veci promjer mlaznice na ulazu
D_ulaz_mlaznice = vanjski_promjer - 2 * debljina_vanjskog_zida_ulaz*1000

debljinaStjenkeRebraNaGrluMlaznice_ulaz_mlaznice = 7 #mm
# fi_alu_rad_ulaz_mlaznice = 2 * debljinaStjenkeRebraNaGrluMlaznice_ulaz_mlaznice / D_ulaz_mlaznice
# fi_alu_stupnjevi_ulaz_mlaznice = np.degrees(fi_alu_rad_ulaz_mlaznice)
fi_alu_rad_ulaz_mlaznice = fi_alu_rad
fi_alu_stupnjevi_ulaz_mlaznice = np.degrees(fi_alu_rad_ulaz_mlaznice)
# print(f"fi_alu_stupnjevi_ulaz_mlaznice: {fi_alu_stupnjevi_ulaz_mlaznice:.6f} ")

# manji promjer mlaznice na ulazu
d_ulaz_mlaznice = np.sqrt(((D_ulaz_mlaznice/1000)**2 * (1 - n *fi_alu_rad_ulaz_mlaznice/(np.pi)) - Auk * (1+k) * 4/np.pi)/(1 - n *fi_alu_rad_ulaz_mlaznice/(np.pi)))
print(f"d_ulaz_mlaznice: {d_ulaz_mlaznice*1000:.6f} mm")

cirkularna_debljina_kanala_ulaza_mlaznice = (D_ulaz_mlaznice - d_ulaz_mlaznice*1000)/2
print(f"cirkularna_debljina_kanala_ulaza_mlaznice: {cirkularna_debljina_kanala_ulaza_mlaznice:.4f} mm\n")

fi_IPA_ulaz_gore_rad = Auk * (8/n) * (1/((D_ulaz_mlaznice/1000)**2 - (d_ulaz_mlaznice)**2))
fi_IPA_ulaz_gore_deg = np.degrees(fi_IPA_ulaz_gore_rad)
duzina_IPA_ulaz_gore_kanala = fi_IPA_ulaz_gore_rad * d_ulaz_mlaznice / 2

print(f"fi_IPA_ulaz_gore_deg: {fi_IPA_ulaz_gore_deg:.6f}°")
print(f"duzina_IPA_ulaz_gore_kanala: {duzina_IPA_ulaz_gore_kanala*1000:.4f} mm")

fi_IPA_izlaz_gore_rad = k * Auk * (8/n) * (1/((D_ulaz_mlaznice/1000)**2 - (d_ulaz_mlaznice)**2))
fi_IPA_izlaz_gore_deg = np.degrees(fi_IPA_izlaz_gore_rad)
duzina_IPA_izlaz_gore_kanala = fi_IPA_izlaz_gore_rad * d_ulaz_mlaznice / 2

print(f"fi_IPA_izlaz_gore_deg: {fi_IPA_izlaz_gore_deg:.6f}°")
print(f"duzina_IPA_izlaz_gore_kanala: {duzina_IPA_izlaz_gore_kanala*1000:.4f} mm\n")

test360 = n *(2 * fi_alu_stupnjevi_ulaz_mlaznice + fi_IPA_ulaz_gore_deg + fi_IPA_izlaz_gore_deg)
print(f"test-ulaz mlaznice-360: {test360:.6f}°\n")


debljina_vanjskog_zida_ulaz = round((faktor_sigurnosti * tlak_fluida* 1e5 * vanjski_promjer / 1e3) / (2*cvrstoca_materijala * 1e6),4) #m
print(f"debljina_vanjskog_zida_ulaz: {debljina_vanjskog_zida_ulaz*1e3} mm")
debljina_unutarnjeg_zida_ulaz = round((faktor_sigurnosti * (tlak_fluida - tlak_komore)* 1e5 * d_ulaz_mlaznice) / (2*cvrstoca_materijala * 1e6),4) #m
print(f"debljina_unutarnjeg_zida_ulaz: {debljina_unutarnjeg_zida_ulaz*1e3} mm\n")

deg_grlo = fi_alu_stupnjevi + fi_IPA_ulaz_deg + fi_alu_stupnjevi + fi_IPA_izlaz_gore_deg
print(f"revolve za deg_grlo: {deg_grlo:.2f} stupnjevi")
print(f"sirina stjenke u stupnjevima za ulaz mlaznice : {fi_alu_stupnjevi_ulaz_mlaznice:.2f} stupnjevi\n")


vanjski_promjer_grlo = d_vani_grlo * 1000 + 2 * debljina_vanjskog_zida_grlo * 1000
unutarnji_promjer_komora = d_ulaz_mlaznice*1000 - 2 * debljina_unutarnjeg_zida_ulaz * 1000
visina_prolaz_izmedju_kanala = 1.5 * (povrsina_cijevi / (n*(d_vani_grlo - unutarnja_stjenka_rashladnih_kanala/1000)))


#dizajn komore
print("dizajn komore izgaranja")
print(f"vanjski promjer: {vanjski_promjer} mm")
print(f"vanjsko-unutarnji promjer - vanjski za rashladne kanale: {D_ulaz_mlaznice:.2f} mm")
print(f"unutarnje-vanjski promjer - unutarnji za rashladne kanale: {d_ulaz_mlaznice*1000:.2f} mm")
print(f"unutarnji_promjer_komora: {unutarnji_promjer_komora:.2f} mm\n")

print("dizajn grla")
print(f"'unutarnji_promjer_grla_mlaznice: {promjer_grla_mlaznice} mm")
print(f"unutarnje-vanjski promjer - unutarnji za rashladne kanale: {unutarnja_stjenka_rashladnih_kanala:.2f} mm")
print(f"vanjsko-unutarnji promjer - vanjski za rashladne kanale: {d_vani_grlo * 1000 :.2f} mm")
print(f"vanjski_promjer_grlo: {vanjski_promjer_grlo:.2f} mm\n")

print(f"visina za _prolaz_izmedju_kanala: {visina_prolaz_izmedju_kanala * 1000:.2f} mm")




# freecad FreeCADCmd kote-mlaznica.py #gui
# freecad --console kote-mlaznica.py  #no gui

import FreeCAD

# 1. Open document
script_dir = os.path.dirname(os.path.realpath(__file__))
# Construct relative path to your document
doc_path = os.path.join(script_dir, "mlaznicaV3-CFD.FCStd")
doc_path = str(doc_path)  # Ensure doc_path is a string
print(f"Opening document: {doc_path}")
doc = FreeCAD.openDocument(doc_path)

#doc = FreeCAD.openDocument("/home/rtx3060/Desktop/Freecad/mlaznicaV3-CFD.FCStd")

# 2. Access spreadsheet
sheet = doc.getObject("Spreadsheet")
if not sheet:
    raise Exception("Spreadsheet 'Params' not found")
    
# 3. Update values using aliases [1][4]
sheet.set('A2', str(vanjski_promjer))
sheet.set('A3', str(D_ulaz_mlaznice))
sheet.set('A4', str(d_ulaz_mlaznice*1000))
sheet.set('A5', str(unutarnji_promjer_komora))

sheet.set('A9', str(promjer_grla_mlaznice))
sheet.set('A10', str(unutarnja_stjenka_rashladnih_kanala))
sheet.set('A11', str(d_vani_grlo * 1000))
sheet.set('A12', str(vanjski_promjer_grlo))

sheet.set('C14', str(deg_grlo))
sheet.set('C15', str(fi_alu_stupnjevi_ulaz_mlaznice))
sheet.set('C16', str(fi_IPA_izlaz_gore_deg))
sheet.set('C17', str(fi_IPA_ulaz_gore_deg))

sheet.set('C19', str(visina_prolaz_izmedju_kanala * 1000))


# 4. Force update chain
sheet.recompute()
doc.recompute()
# FreeCADGui.updateGui()  # Refresh if GUI exists [2]


# 6. Automatsko generiranje STL za FaceBindere
import os
import Mesh
from FreeCAD import Base
import re  

def export_facebinders_to_stl():
    output_dir = os.path.join(os.path.dirname(doc.FileName), "cases")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        # Remove existing .stl files in the cases folder
    for file_name in os.listdir(output_dir):
        print(file_name)
        if file_name.endswith(".stl"):
            file_path = os.path.join(output_dir, file_name)
            try:
                os.remove(file_path)
                print(f"Deleted: {file_path}")
            except Exception as e:
                print(f"Failed to delete {file_path}: {str(e)}")
    success_count = 0
    error_count = 0
    total_facets = 0

    # Memory optimization
    import gc
    import subprocess
    import shutil
    gc.collect()

    for obj in doc.Objects:
        try:
            if "facebinder" not in obj.Name.lower():
                continue

            sanitized_label = re.sub(r'[^\w-]', '_', obj.Label).strip('_') or obj.Name
            print(f"\nProcessing: {sanitized_label}")

            if not obj.Shape or obj.Shape.isNull():
                print("  - Skipped: Empty shape")
                error_count += 1
                continue

            # Adjusted tessellation tolerance
            facets = obj.Shape.tessellate(1e-2)  # Reduced precision for memory
            if not facets[1]:
                print("  - Skipped: No facets")
                error_count += 1
                continue

            stl_path = os.path.join(output_dir, f"{sanitized_label}.stl")
            with open(stl_path, 'w') as f:
                f.write(f"solid {sanitized_label}\n")
                
                # Facet processing with safety checks
                facet_count = 0
                for triangle in facets[1]:
                    try:
                        v1, v2, v3 = [facets[0][i] for i in triangle]
                        
                        # Normal calculation with validation
                        normal = (v2 - v1).cross(v3 - v1)
                        if normal.Length < 1e-10:
                            continue  # Skip degenerate facets
                            
                        normal.normalize()
                        
                        f.write(
                            f"  facet normal {normal.x:.7e} {normal.y:.7e} {normal.z:.7e}\n"
                            "    outer loop\n"
                            f"      vertex {v1.x:.7e} {v1.y:.7e} {v1.z:.7e}\n"
                            f"      vertex {v2.x:.7e} {v2.y:.7e} {v2.z:.7e}\n"
                            f"      vertex {v3.x:.7e} {v3.y:.7e} {v3.z:.7e}\n"
                            "    endloop\n"
                            "  endfacet\n"
                        )
                        facet_count += 1
                    except Exception as e:
                        print(f"  - Skipped facet: {str(e)}")

                f.write(f"endsolid {sanitized_label}\n")
                total_facets += facet_count
                print(f"  - Exported {facet_count} triangles")

            success_count += 1
            gc.collect()  # Manual garbage collection

        except Exception as e:
            error_count += 1
            print(f"ERROR: {str(e)}")
            FreeCAD.Console.PrintError(f"{sanitized_label} error: {str(e)}\n")

    print(f"\nExport complete: {success_count} files, {total_facets} triangles total")
    # Execute shell commands to concatenate STL files

    output_dir = os.path.join(os.path.dirname(doc.FileName), "cases")
    solid_domena_path = os.path.join(output_dir, "solid-domena.stl")
    fluid_domena_path = os.path.join(output_dir, "fluid-domena.stl")

    # Concatenate files for solid-domena.stl
    # Verify existence of required STL files before concatenating
    required_solid_files = [
        "intersection-wall.stl", "wall-unutra.stl", "wall-unutra-ref.stl",
        "wall-vani.stl", "wall-vani-ref.stl", "wedge-wall-coolant.stl", "wedge-wall.stl"
    ]
    required_fluid_files = [
        "inlet.stl", "outlet.stl", "intersection-wall.stl", "wall-outlet.stl"
    ]

    def concatenate_files(file_list, output_path):
        with open(output_path, 'w') as outfile:
            for file_name in file_list:
                file_path = os.path.join(output_dir, file_name)
                if os.path.exists(file_path):
                    with open(file_path, 'r') as infile:
                        outfile.write(infile.read())
                else:
                    print(f"WARNING: {file_path} does not exist and will be skipped.")

    # Concatenate files for solid-domena.stl
    concatenate_files(required_solid_files, solid_domena_path)

    # Concatenate files for fluid-domena.stl
    concatenate_files(required_fluid_files, fluid_domena_path)

    # Delete old .stl files in cases/solid-domena/ and cases/fluid-domena/
    solid_domena_dir = os.path.join(output_dir, "solid-domena")
    fluid_domena_dir = os.path.join(output_dir, "fluid-domena")

    for folder in [solid_domena_dir, fluid_domena_dir]:
        if os.path.exists(folder):
            for file_name in os.listdir(folder):
                if file_name.endswith(".stl"):
                    file_path = os.path.join(folder, file_name)
                    try:
                        os.remove(file_path)
                        print(f"Deleted: {file_path}")
                    except Exception as e:
                        print(f"Failed to delete {file_path}: {str(e)}")

    # Copy the new solid-domena.stl and fluid-domena.stl to their respective directories
    os.makedirs(solid_domena_dir, exist_ok=True)
    os.makedirs(fluid_domena_dir, exist_ok=True)

    solid_domena_target = os.path.join(solid_domena_dir, "solid-domena.stl")
    fluid_domena_target = os.path.join(fluid_domena_dir, "fluid-domena.stl")

    try:
        shutil.copy(solid_domena_path, solid_domena_target)
        print(f"Copied {solid_domena_path} to {solid_domena_target}")
        shutil.copy(fluid_domena_path, fluid_domena_target)
        print(f"Copied {fluid_domena_path} to {fluid_domena_target}")
    except Exception as e:
        print(f"Failed to copy files: {str(e)}")

    return success_count > 0




# Run with memory protection
try:
    export_facebinders_to_stl()
except MemoryError:
    print("\nFATAL: Out of memory - try reducing model complexity")
except Exception as e:
    print(f"\nCritical error: {str(e)}")
finally:
    exit()

# cat intersection-wall.stl  wall-unutra.stl wall-unutra-ref.stl wall-vani.stl wall-vani-ref.stl wedge-wall-coolant.stl wedge-wall.stl >> solid-domena.stl
# cat inlet.stl outlet.stl intersection-wall.stl wall-outlet.stl >> fluid-domena.stl    