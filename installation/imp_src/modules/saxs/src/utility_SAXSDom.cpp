/**
 *  \file IMP/saxs/utility.cpp
 *  \brief Functions to deal with very common saxs operations
 *  Copyright 2007-2016 IMP Inventors. All rights reserved.
*/

//#include <IMP/saxs/utility.h>
#include <IMP/saxs/utility_SAXSDom.h>
#include <IMP/saxs/SolventAccessibleSurface.h>
#include <IMP/atom/pdb_SAXSDom.h>

IMPSAXS_BEGIN_NAMESPACE

/*
Profile* compute_profile(Particles particles, double min_q,
                         double max_q, double delta_q, FormFactorTable* ft,
                         FormFactorType ff_type, bool hydration_layer, bool fit,
                         bool reciprocal, bool ab_initio, bool vacuum,
                         std::string beam_profile_file) {
  IMP_NEW(Profile, profile, (min_q, max_q, delta_q));
  if (reciprocal) profile->set_ff_table(ft);
  if (beam_profile_file.size() > 0) profile->set_beam_profile(beam_profile_file);

  // compute surface accessibility and average radius
  Vector<double> surface_area;
  SolventAccessibleSurface s;
  double average_radius = 0.0;
  if (hydration_layer) {
    // add radius
    for (unsigned int i = 0; i < particles.size(); i++) {
      double radius = ft->get_radius(particles[i], ff_type);
      core::XYZR::setup_particle(particles[i], radius);
      average_radius += radius;
    }
    surface_area = s.get_solvent_accessibility(core::XYZRs(particles));
    average_radius /= particles.size();
    profile->set_average_radius(average_radius);
  }

  // pick profile calculation based on input parameters
  if (!fit) {         // regular profile, no c1/c2 fitting
    if (ab_initio) {  // bead model, constant form factor
      profile->calculate_profile_constant_form_factor(particles);
    } else if (vacuum) {
      profile->calculate_profile_partial(particles, surface_area, ff_type);
      profile->sum_partial_profiles(0.0, 0.0);  // c1 = 0;
    } else {
      profile->calculate_profile(particles, ff_type, reciprocal);
    }
  } else {  // c1/c2 fitting
    if (reciprocal)
      profile->calculate_profile_reciprocal_partial(particles, surface_area,
                                                    ff_type);
    else
      profile->calculate_profile_partial(particles, surface_area, ff_type);
  }
  return profile.release();
}
*/

/*
void read_pdb(const std::string file, std::vector<std::string>& pdb_file_names,
              std::vector<IMP::Particles>& particles_vec,
              bool residue_level, bool heavy_atoms_only, int multi_model_pdb) {

  IMP::Model* model = new IMP::Model();

  IMP::atom::Hierarchies mhds;
  IMP::atom::PDBSelector* selector;
  if (residue_level)  // read CA only
    selector = new IMP::atom::CAlphaPDBSelector();
  else if (heavy_atoms_only)  // read without hydrogens
    selector = new IMP::atom::NonWaterNonHydrogenPDBSelector();
  else  // read with hydrogens
    selector = new IMP::atom::NonWaterPDBSelector();

  if (multi_model_pdb == 2) {
    mhds = read_multimodel_pdb(file, model, selector, true);
  } else {
    if (multi_model_pdb == 3) {
      IMP::atom::Hierarchy mhd =
          IMP::atom::read_pdb(file, model, selector, false, true);
      mhds.push_back(mhd);
    } else {
      IMP::atom::Hierarchy mhd =
          IMP::atom::read_pdb(file, model, selector, true, true);
      mhds.push_back(mhd);
    }
  }

  for (unsigned int h_index = 0; h_index < mhds.size(); h_index++) {
    IMP::ParticlesTemp ps =
        get_by_type(mhds[h_index], IMP::atom::ATOM_TYPE);
    if (ps.size() > 0) {  // pdb file
      std::string pdb_id = file;
      if (mhds.size() > 1) {
        pdb_id = trim_extension(file) + "_m" +
                 std::string(boost::lexical_cast<std::string>(h_index + 1)) +
                 ".pdb";
      }
      pdb_file_names.push_back(pdb_id);
      particles_vec.push_back(IMP::get_as<IMP::Particles>(ps));
      std::cout << ps.size() << " atoms were read from PDB file " << file;
      if (mhds.size() > 1) std::cout << " MODEL " << h_index + 1;
      std::cout << std::endl;
    }
  }
}
*/

// added by jie 11/21/2016
void read_pdb_saxs(const std::string file, std::vector<std::string>& pdb_file_names, std::string inputString,
              std::vector<IMP::Particles>& particles_vec,
              bool residue_level, bool heavy_atoms_only, int multi_model_pdb) {

  IMP::Model* model = new IMP::Model();

  IMP::atom::Hierarchies mhds;
  IMP::atom::PDBSelector* selector;
  if (residue_level)  // read CA only
    selector = new IMP::atom::CAlphaPDBSelector();
  else if (heavy_atoms_only)  // read without hydrogens
    selector = new IMP::atom::NonWaterNonHydrogenPDBSelector();
  else  // read with hydrogens
    selector = new IMP::atom::NonWaterPDBSelector();

  if (multi_model_pdb == 2) {
    mhds = read_multimodel_pdb(file, model, selector, true);
  } else {
    if (multi_model_pdb == 3) {
      IMP::atom::Hierarchy mhd =
          IMP::atom::read_pdb_saxs(file, model,inputString, selector, false, true);
      mhds.push_back(mhd);
    } else {
      IMP::atom::Hierarchy mhd =
          IMP::atom::read_pdb_saxs(file, model,inputString, selector, true, true);
      mhds.push_back(mhd);
    }
  }
  for (unsigned int h_index = 0; h_index < mhds.size(); h_index++) {
    IMP::ParticlesTemp ps =
        get_by_type(mhds[h_index], IMP::atom::ATOM_TYPE);
    if (ps.size() > 0) {  // pdb file
      std::string pdb_id = file;
      if (mhds.size() > 1) {
        pdb_id = trim_extension(file) + "_m" +
                 std::string(boost::lexical_cast<std::string>(h_index + 1)) +
                 ".pdb";
      }
      pdb_file_names.push_back(pdb_id);
      particles_vec.push_back(IMP::get_as<IMP::Particles>(ps));
      //std::cout << ps.size() << " atoms were read from PDB file " << file;
      if (mhds.size() > 1) std::cout << " MODEL " << h_index + 1;
    }
  }
}

/*
void read_files(const std::vector<std::string>& files,
                std::vector<std::string>& pdb_file_names,
                std::vector<std::string>& dat_files,
                std::vector<IMP::Particles>& particles_vec,
                Profiles& exp_profiles, bool residue_level,
                bool heavy_atoms_only, int multi_model_pdb, float max_q) {

  for (unsigned int i = 0; i < files.size(); i++) {
    // check if file exists
    std::ifstream in_file(files[i].c_str());
    if (!in_file) {
      std::cerr << "Can't open file " << files[i] << std::endl;
      return;
    }
    // 1. try as pdb
    try {
      read_pdb(files[i], pdb_file_names, particles_vec, residue_level,
               heavy_atoms_only, multi_model_pdb);
    }
    catch (IMP::ValueException e) {  // not a pdb file
      // 2. try as a dat profile file
      IMP_NEW(Profile, profile, (files[i], false, max_q));
      if (profile->size() == 0) {
        std::cerr << "can't parse input file " << files[i] << std::endl;
        return;
      } else {
        dat_files.push_back(files[i]);
        exp_profiles.push_back(profile);
        std::cout << "Profile read from file " << files[i]
                  << " size = " << profile->size() << std::endl;
      }
    }
  }
}
*/

void read_files_saxs(const std::vector<std::string>& files,
                std::vector<std::string>& pdb_file_names,
                std::vector<std::string>& dat_files, std::string inputString,
                std::vector<IMP::Particles>& particles_vec,
                Profiles& exp_profiles, bool residue_level,
                bool heavy_atoms_only, int multi_model_pdb, float max_q) {
   
  std::vector<std::string> files_pdb, files_dat;
  std::string pdb_subfix (".pdb");
  std::string dat_subfix (".dat");
  for (unsigned int i = 0; i < files.size(); i++) {
	  if (files[i].find(pdb_subfix) != std::string::npos) 
	  {
		files_pdb.push_back(files[i]);
	  }else if(files[i].find(dat_subfix) != std::string::npos) 
	  {
		  files_dat.push_back(files[i]);
	  }else{
		  std::cout << "Wrong format files: "<< files[i] << std::endl;
		  exit(-1);
	  }
  }
  if(files_pdb.size() > 1 or files_dat.size() >1)
  {
	  std::cout << "Currently one 1 pdb and 1 dat is accepted"<< std::endl;
	  exit(-1);
  }
  /*
  for (unsigned int i = 0; i < files.size(); i++) {
    // check if file exists
    //std::ifstream in_file(files[i].c_str());
    //if (!in_file) {
    //  std::cerr << "Can't open file " << files[i] << std::endl;
    //  return;
    //}
    // 1. try as pdb
    try {
	  //std::cout << "Try pdb in utility_saxs_jie.cpp " << files[i] << std::endl;
      read_pdb_saxs(files[i], pdb_file_names,inputString, particles_vec, residue_level,
               heavy_atoms_only, multi_model_pdb);
    }
    catch (IMP::ValueException e) {  // not a pdb file
      // 2. try as a dat profile file
	  //std::cout << "Try dat in utility_saxs_jie.cpp " << files[i] << std::endl;
      IMP_NEW(Profile, profile, (files[i], false, max_q));
      if (profile->size() == 0) {
        std::cerr << "can't parse input file " << files[i] << std::endl;
        return;
      } else {
        dat_files.push_back(files[i]);
        exp_profiles.push_back(profile);
        std::cout << "Profile read from file " << files[i]
                  << " size = " << profile->size() << std::endl;
      }
    }
  }
  */
  for (unsigned int i = 0; i < files_pdb.size(); i++) {
	  //std::cout << "Try pdb in utility_saxs_jie.cpp " << files[i] << std::endl;
      read_pdb_saxs(files_pdb[i], pdb_file_names,inputString, particles_vec, residue_level,
               heavy_atoms_only, multi_model_pdb);
    
  }
  for (unsigned int i = 0; i < files_dat.size(); i++) {
      // 2. try as a dat profile file
	  //std::cout << "Try dat in utility_saxs_jie.cpp " << files[i] << std::endl;
      IMP_NEW(Profile, profile, (files_dat[i], false, max_q));
      if (profile->size() == 0) {
        std::cerr << "can't parse input file " << files_dat[i] << std::endl;
        return;
      } else {
        dat_files.push_back(files_dat[i]);
        exp_profiles.push_back(profile);
        //std::cout << "Profile read from file " << files_dat[i]
                //  << " size = " << profile->size() << std::endl;
      }
  }
}



void read_files_saxs_pdbOnly(const std::vector<std::string>& files,
                std::vector<std::string>& pdb_file_names,
                std::vector<std::string>& dat_files, std::string inputString,
                std::vector<IMP::Particles>& particles_vec,
                bool residue_level,
                bool heavy_atoms_only, int multi_model_pdb, float max_q) {
   
  std::vector<std::string> files_pdb, files_dat;
  std::string pdb_subfix (".pdb");
  std::string dat_subfix (".dat");
  for (unsigned int i = 0; i < files.size(); i++) {
	  if (files[i].find(pdb_subfix) != std::string::npos) 
	  {
		files_pdb.push_back(files[i]);
	  }else if(files[i].find(dat_subfix) != std::string::npos) 
	  {
		  files_dat.push_back(files[i]);
	  }else{
		  std::cout << "Wrong format files: "<< files[i] << std::endl;
		  exit(-1);
	  }
  }
  if(files_pdb.size() > 1 or files_dat.size() >1)
  {
	  std::cout << "Currently one 1 pdb and 1 dat is accepted"<< std::endl;
	  exit(-1);
  }
  /*
  for (unsigned int i = 0; i < files.size(); i++) {
    // check if file exists
    //std::ifstream in_file(files[i].c_str());
    //if (!in_file) {
    //  std::cerr << "Can't open file " << files[i] << std::endl;
    //  return;
    //}
    // 1. try as pdb
    try {
	  //std::cout << "Try pdb in utility_saxs_jie.cpp " << files[i] << std::endl;
      read_pdb_saxs(files[i], pdb_file_names,inputString, particles_vec, residue_level,
               heavy_atoms_only, multi_model_pdb);
    }
    catch (IMP::ValueException e) {  // not a pdb file
      // 2. try as a dat profile file
	  //std::cout << "Try dat in utility_saxs_jie.cpp " << files[i] << std::endl;
      IMP_NEW(Profile, profile, (files[i], false, max_q));
      if (profile->size() == 0) {
        std::cerr << "can't parse input file " << files[i] << std::endl;
        return;
      } else {
        dat_files.push_back(files[i]);
        exp_profiles.push_back(profile);
        std::cout << "Profile read from file " << files[i]
                  << " size = " << profile->size() << std::endl;
      }
    }
  }
  */
  for (unsigned int i = 0; i < files_pdb.size(); i++) {
	  //std::cout << "Try pdb in utility_saxs_jie.cpp " << files[i] << std::endl;
      read_pdb_saxs(files_pdb[i], pdb_file_names,inputString, particles_vec, residue_level,
               heavy_atoms_only, multi_model_pdb);
    
  }
}


/*
std::string trim_extension(const std::string file_name) {
  if (file_name[file_name.size() - 4] == '.')
    return file_name.substr(0, file_name.size() - 4);
  return file_name;
}
*/
IMPSAXS_END_NAMESPACE
