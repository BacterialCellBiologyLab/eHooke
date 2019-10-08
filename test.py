from ehooke import EHooke

app = EHooke()
app.load_base_image("test_phase.tif")
app.compute_mask()
app.load_fluor_image("test_membrane.tif")
app.load_option_image("test_dna.tif")
app.compute_segments()
app.compute_cells()
app.parameters.cellprocessingparams.find_septum = True
app.parameters.cellprocessingparams.classify_cells = True
app.parameters.cellprocessingparams.microscope = "Epifluorescence"
app.process_cells()
app.generate_reports()
app.compute_coloc()