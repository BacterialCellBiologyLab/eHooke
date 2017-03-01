from ehooke import EHooke

t = EHooke()
t.load_base_image("test_phase.tif")
t.compute_mask()
t.load_fluor_image("test_fluor.tif")
t.compute_segments()
t.compute_cells()
t.merge_cells(340, 346)
t.split_cells(346)
t.process_cells()
