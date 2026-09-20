RSV-selection-pressure
Analysis code for selection-pressure inference and coordinate reconciliation in the RSV-A and RSV-B F and G genes, comparing pre- and post-rollout periods.
Software
HyPhy v2.5 (FEL, FUBAR, MEME), MAPLE v0.7.5, Nextclade, PyMOL, Python 3. Analyses were run on Linux (Ubuntu under WSL2).
Scripts
selection_pressure_pipeline_6.py — gene extraction with empirical reading-frame determination, quality control, collapsing to unique haplotypes with matched tree pruning, HyPhy orchestration, and parsing of results. Run with --status, --extract, --run, --parse.
find_offset.py — determines the offset between analysed alignment columns and reference F0 numbering by global alignment of the consensus translation against the subtype-specific reference.
map_G_positions.py — equivalent mapping for the G gene, which requires alignment rather than a single offset because of the ON1 and BA genotype duplications.
rsv_F_domain_map.py / rsv_G_domain_map.py — linear domain maps of positively selected sites.
rsv_F_surface.py — structural mapping onto prefusion F trimers (PDB 4MMV for RSV-A, 5UDE for RSV-B) in PyMOL.
Coordinate systems
Positions returned by the selection analysis are alignment column indices, not reference protein positions. Offsets of +2 (RSV-A F and G), +6 (RSV-B F) and +5 (RSV-B G) convert them to reference numbering.
Data
No sequence data is included. Sequences were obtained from GISAID under EPI_SET_260918mg (https://doi.org/10.55876/gis8.260918mg).
License
MIT
