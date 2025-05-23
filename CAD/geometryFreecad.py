import os
import json
import FreeCAD

# Load the design data from the JSON file
script_dir = os.path.dirname(os.path.realpath(__file__))
data_file = os.path.join(script_dir, "design_data.json")
with open(data_file, "r") as f:
    regenNozzleDesign = json.load(f)

print(f"Loaded design data: {regenNozzleDesign}")

# FreeCAD operations
# doc_path = os.path.join(script_dir, "mlaznica3DPrint.FCStd")
doc_path = os.path.join(script_dir, regen_data["fcstd_filename"])
doc = FreeCAD.openDocument(doc_path)

sheet = doc.getObject("Spreadsheet")
if not sheet:
    raise Exception("Spreadsheet 'Params' not found")

# Update spreadsheet values
sheet.set('A2', str(regenNozzleDesign["D_inlet"] * 1000))
sheet.set('A3', str(regenNozzleDesign["D_ulaz_mlaznice"] * 1000))
sheet.set('A4', str(regenNozzleDesign["d_cooling_channels_inlet"] * 1000))
sheet.set('A5', str(regenNozzleDesign["combustion_chamber_diameter"] * 1000))

sheet.set('A9', str(regenNozzleDesign["promjer_grla_mlaznice"]))
sheet.set('A10', str(regenNozzleDesign["d_cooling_channels_throat"] * 1000))
sheet.set('A11', str(regenNozzleDesign["D_cooling_channels_throat"] * 1000))
sheet.set('A12', str(regenNozzleDesign["D_throat"] * 1000))

sheet.set('C14', str(regenNozzleDesign["fi_revolve"]))
sheet.set('C15', str(regenNozzleDesign["fi_alu_stupnjevi"]))
sheet.set('C16', str(regenNozzleDesign["fi_inlet_stupnjevi_ulaz"]))
sheet.set('C17', str(regenNozzleDesign["fi_inlet_stupnjevi_izlaz"]))
sheet.set('C19', str(regenNozzleDesign["channel_height"] * 1000))

sheet.set('C21', str(regenNozzleDesign["n_revolve"]))

sheet.recompute()
doc.recompute()