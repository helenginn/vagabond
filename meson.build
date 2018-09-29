project('vagabond', 'cpp', 'c')
qt5 = import('qt5')

qt5_dep = dependency('qt5', modules: ['Core', 'Gui', 'Widgets'], required : false)
dep_gl = dependency('gl', required : false)

png_dep = dependency('libpng')
boost_dep = dependency('boost')
fftw_dep = dependency('fftw3f')

# Geometry and/or rotamer lookup madness

lib_vgeom = library('vgeom', 'libinfo/GeomTable.cpp', 'libinfo/RotamerTable.cpp', cpp_args: ['-w'])
inc_vgeom = include_directories('libinfo')
dep_vgeom = declare_dependency(link_with: lib_vgeom, include_directories: inc_vgeom)

# numerical recipes, single value decomposition

lib_ica = library('ica', 'libica/matrix.cpp', 'libica/svdcmp.cpp', 'liblbfgs/lbfgs.c', cpp_args: ['-w'])
inc_ica = include_directories('libica')
dep_ica = declare_dependency(link_with: lib_ica, include_directories: inc_ica)

# cmtzlib / csymlib: SYMINFO-variable-independent distributions needs separate compiling
# May or may not need the maths library on various systems

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)

lib_ccp4 = library('ccp4', 'libccp4/ccp4_array.c', 'libccp4/ccp4_general.c', 'libccp4/ccp4_parser.c', 'libccp4/ccp4_program.c', 'libccp4/ccp4_unitcell.c', 'libccp4/cmap_accessor.c', 'libccp4/cmap_close.c', 'libccp4/cmap_data.c', 'libccp4/cmap_header.c', 'libccp4/cmap_labels.c', 'libccp4/cmap_open.c', 'libccp4/cmap_skew.c', 'libccp4/cmap_stats.c', 'libccp4/cmap_symop.c', 'libccp4/cmtzlib.c', 'libccp4/csymlib.c', 'libccp4/cvecmat.c', 'libccp4/library_err.c', 'libccp4/library_file.c', 'libccp4/library_utils.c', 'libccp4/pack_c.c', dependencies: m_dep, c_args: ['-w'])
inc_ccp4 = include_directories('libccp4')
dep_ccp4 = declare_dependency(link_with: lib_ccp4, include_directories: inc_ccp4)

# Meat of the vagabond library, minus command line-specific or GUI-specific source files.

lib_vag= library('vag', 'libsrc/Absolute.cpp', 'libsrc/Anchor.cpp', 'libsrc/Anisotropicator.cpp', 'libsrc/Atom.cpp', 'libsrc/AtomGroup.cpp', 'libsrc/Backbone.cpp', 'libsrc/Bond.cpp', 'libsrc/BondGroup.cpp', 'libsrc/BucketBulkSolvent.cpp', 'libsrc/Bucket.cpp', 'libsrc/BucketMonteCarlo.cpp', 'libsrc/Chromosomal.cpp', 'libsrc/Crystal.cpp', 'libsrc/CSV.cpp', 'libsrc/Diffraction.cpp', 'libsrc/DiffractionMTZ.cpp', 'libsrc/Distributor.cpp', 'libsrc/Element.cpp', 'libsrc/ExplicitModel.cpp', 'libsrc/fftw3d.cpp', 'libsrc/FileReader.cpp', 'libsrc/FlexGlobal.cpp', 'libsrc/FlexLocal.cpp', 'libsrc/Hydrogenator.cpp', 'libsrc/Knotter.cpp', 'libsrc/LeftThing.cpp', 'libsrc/mat3x3.cpp', 'libsrc/mat4x4.cpp', 'libsrc/maths.cpp', 'libsrc/MDNode.cpp', 'libsrc/Model.cpp', 'libsrc/Molecule.cpp', 'libsrc/Monomer.cpp', 'libsrc/Options.cpp', 'libsrc/ParamBand.cpp', 'libsrc/Parser.cpp', 'libsrc/PDBReader.cpp', 'libsrc/Plucker.cpp', 'libsrc/PNGFile.cpp', 'libsrc/Polymer.cpp', 'libsrc/RefinementGridSearch.cpp', 'libsrc/RefinementLBFGS.cpp', 'libsrc/RefinementNelderMead.cpp', 'libsrc/RefinementStepSearch.cpp', 'libsrc/RefinementStrategy.cpp', 'libsrc/Reflex.cpp', 'libsrc/Sampler.cpp', 'libsrc/Shouter.cpp', 'libsrc/Sidechain.cpp', 'libsrc/SSRigger.cpp', 'libsrc/StateValue.cpp', 'libsrc/TextManager.cpp', 'libsrc/Thing.cpp', 'libsrc/Timer.cpp', 'libsrc/VBondReader.cpp', 'libsrc/vec3.cpp', 'libsrc/VScope.cpp', 'libsrc/VScript.cpp', 'libsrc/WaterCluster.cpp', 'libsrc/WaterNetwork.cpp', dependencies: [png_dep, fftw_dep, dep_vgeom, dep_ica, dep_ccp4, boost_dep], cpp_args: ['-w', '-std=c++11'])
inc_vag = include_directories('libsrc')
dep_vag = declare_dependency(link_with: lib_vag, include_directories: inc_vag)

# Command line interface

executable('vagabond', 'libsrc/main.cpp', dependencies: [dep_vag, boost_dep], install: true)

# Graphical user interface
# (including all the Qt and OpenGL stuff)

# Pre-processing by Qt5

if qt5_dep.found() and dep_gl.found()
	moc_files = qt5.preprocess(moc_headers : ['libgui/qtgui/CrystalExplorer.h', 'libgui/qtgui/Dialogue.h', 'libgui/qtgui/MoleculeExplorer.h', 'libgui/qtgui/MonomerExplorer.h', 'libgui/qtgui/ResButton.h', 'libgui/qtgui/SequenceView.h', 'libgui/qtgui/SetterEdit.h', 'libgui/qtgui/StartScreen.h', 'libgui/qtgui/VagWindow.h'],
			moc_extra_arguments: ['-DMAKES_MY_MOC_HEADER_COMPILE'])

# Non-Qt5 GUI elements (mostly OpenGL stuff)
# Uses Qt5-supplied OpenGL

	lib_vgui = library('vgui', 'libgui/Density2GL.cpp', 'libgui/GLKeeper.cpp', 'libgui/GLObject.cpp', 'libgui/shader.cpp', 'libgui/Vagabond2GL.cpp', cpp_args: ['-w'], dependencies : [ qt5_dep, dep_vag ])
	inc_vgui = include_directories('libgui')
dep_vgui = declare_dependency(link_with: lib_vgui, include_directories: inc_vgui)

	executable('vagabond-gui', 'libgui/qtgui/gui_main.cpp', 'libgui/qtgui/CrystalExplorer.cpp', 'libgui/qtgui/Dialogue.cpp', 'libgui/qtgui/gui_main.cpp', 'libgui/qtgui/InstructionThread.cpp', 'libgui/qtgui/MoleculeExplorer.cpp', 'libgui/qtgui/MonomerExplorer.cpp', 'libgui/qtgui/ResButton.cpp', 'libgui/qtgui/SequenceView.cpp', 'libgui/qtgui/SetterEdit.cpp', 'libgui/qtgui/StartScreen.cpp', 'libgui/qtgui/VagabondGLWidget.cpp', 'libgui/qtgui/VagWindow.cpp', moc_files, dependencies: [qt5_dep, dep_vag, dep_vgeom, dep_vgui], install: true, cpp_args: ['-w', '-std=c++11'])
else
	message('QT5 not found, so not building the GUI.')
endif
