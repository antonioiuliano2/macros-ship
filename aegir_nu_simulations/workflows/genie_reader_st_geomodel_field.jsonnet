local lib = import 'lib.libsonnet';
{
  driver: lib.driver(1000),
  sources: {
    field: lib.field_placements(),
    geometry: lib.geomodel_geometry {
      db_file: 'ship_geometry.db',
      sensitive_volumes: [{pattern:'neutrino_detector/hcal/tile',detector:'UpstreamTagger'},
                        {pattern: 'neutrino_detector/mgc/Si',detector:'SurroundTagger'},
                        {pattern:'neutrino_detector/target/Si_X',detector:'TimingDetector'}],
    },
    genie_reader: lib.genie_reader{
      file: '2026_11_16_4e18_NeutrinoDetector_allflavours.rootracker.root'
    },
  },
  modules: {
    // Scoring mode: every step that deposits energy in a sensitive volume
    // becomes a SimHit, which is what Shannon digitises. Tracks below the
    // energy cut are killed; secondaries below particle_ke_cut are tracked
    // but not written as SimParticles.
    // Crossing mode: a sensitive volume records the first step of each
    // track that enters it. Tracks below the energy cut are killed.
    geant4: lib.geant4{
      cpp: 'geant4_module',
      physics_list: 'FTFP_BERT',
      sd_mode: 'scoring',
      energy_cut_threshold: 0.001, // GeV (kill tracks below 1 MeV)
      energy_cut: true,
      particle_ke_cut: 0.1 //GeV (the usual do not store if Ekin<100 MeV cut)
    },
    output: lib.full_output('nu_field_genie_reader_st_output.root', 'nu_field_genie_reader_st_validation.root'),
  },
}
