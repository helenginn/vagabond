project('vagabond', 'cpp', 'c')
qt5 = import('qt5')

qt5_dep = dependency('qt5', modules: ['Core', 'Gui', 'Widgets'], required : false)
dep_gl = dependency('gl', required : false)

png_dep = dependency('libpng')
boost_dep = dependency('boost')
fftw_dep = dependency('fftw3f')

# Geometry and/or rotamer lookup madness

lib_vgeom = library('vgeom', 'libinfo/GeomTable.cpp', 'libinfo/RotamerTable.cpp', cpp_args: ['-w'], install: true)
inc_vgeom = include_directories('libinfo')
dep_vgeom = declare_dependency(link_with: lib_vgeom, include_directories: inc_vgeom)

# numerical recipes, single value decomposition

lib_ica = library('ica', 'libica/matrix.cpp', 'libica/svdcmp.cpp', 'liblbfgs/lbfgs.c', cpp_args: ['-w'], install: true)
inc_ica = include_directories('libica')
dep_ica = declare_dependency(link_with: lib_ica, include_directories: inc_ica)

# cmtzlib / csymlib: SYMINFO-variable-independent distributions needs separate compiling
# May or may not need the maths library on various systems

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)

lib_ccp4 = library('ccp4', 'libccp4/ccp4_array.c', 'libccp4/ccp4_general.c', 'libccp4/ccp4_parser.c', 'libccp4/ccp4_program.c', 'libccp4/ccp4_unitcell.c', 'libccp4/cmap_accessor.c', 'libccp4/cmap_close.c', 'libccp4/cmap_data.c', 'libccp4/cmap_header.c', 'libccp4/cmap_labels.c', 'libccp4/cmap_open.c', 'libccp4/cmap_skew.c', 'libccp4/cmap_stats.c', 'libccp4/cmap_symop.c', 'libccp4/cmtzlib.c', 'libccp4/csymlib.c', 'libccp4/cvecmat.c', 'libccp4/library_err.c', 'libccp4/library_file.c', 'libccp4/library_utils.c', 'libccp4/pack_c.c', dependencies: m_dep, c_args: ['-w'], install: true)
inc_ccp4 = include_directories('libccp4')
dep_ccp4 = declare_dependency(link_with: lib_ccp4, include_directories: inc_ccp4)

# Run the version and hash generator
hash_script = find_program('get_hash.sh')

# This will produce the dependency file.
run_command('get_hash.sh')

gen_src = custom_target('commit-output',
                        input : ['libsrc/commit.h.pre'],
                        output : ['commit.h'],
                        command : [hash_script, '@INPUT@', '@OUTPUT@'],
                        build_by_default: true,
                        build_always_stale: true)


# Meat of the vagabond library, minus command line-specific or GUI-specific source files.

base_files = files('libsrc/Absolute.cpp', 
'libsrc/Anchor.cpp', 
'libsrc/Anisotropicator.cpp', 
'libsrc/Atom.cpp', 
'libsrc/AtomGroup.cpp', 
'libsrc/Backbone.cpp', 
'libsrc/Balance.cpp', 
'libsrc/BaseParser.cpp', 
'libsrc/Bond.cpp', 
'libsrc/BondGroup.cpp', 
'libsrc/BucketBulkSolvent.cpp', 
'libsrc/BucketPerStrand.cpp', 
'libsrc/Bucket.cpp', 
'libsrc/Chelate.cpp', 
'libsrc/Converter.cpp', 
'libsrc/Crystal.cpp', 
'libsrc/CSV.cpp', 
'libsrc/Diffraction.cpp', 
'libsrc/DiffractionMTZ.cpp', 
'libsrc/Element.cpp', 
'libsrc/ExplicitModel.cpp', 
'libsrc/FFT.cpp', 
'libsrc/Fibonacci.cpp', 
'libsrc/FileReader.cpp', 
'libsrc/FlexGlobal.cpp', 
'libsrc/FlexLocal.cpp', 
'libsrc/GhostBond.cpp',
'libsrc/Hydrogenator.cpp', 
'libsrc/KeyPoints.cpp', 
'libsrc/Knotter.cpp', 
'libsrc/mat3x3.cpp', 
'libsrc/mat4x4.cpp', 
'libsrc/maths.cpp', 
'libsrc/Model.cpp', 
'libsrc/Molecule.cpp', 
'libsrc/Motion.cpp', 
'libsrc/Monomer.cpp', 
'libsrc/Novalent.cpp', 
'libsrc/Options.cpp', 
'libsrc/ParamBand.cpp', 
'libsrc/Parser.cpp', 
'libsrc/PartialStructure.cpp', 
'libsrc/PDBReader.cpp', 
'libsrc/PNGFile.cpp', 
'libsrc/Polymer.cpp', 
'libsrc/Refitter.cpp', 
'libsrc/RefinementGridSearch.cpp', 
'libsrc/RefinementLBFGS.cpp', 
'libsrc/RefinementList.cpp', 
'libsrc/RefinementNelderMead.cpp', 
'libsrc/RefinementStepSearch.cpp', 
'libsrc/RefinementStrategy.cpp', 
'libsrc/Sampler.cpp', 
'libsrc/Shouter.cpp', 
'libsrc/Sidechain.cpp', 
'libsrc/SpaceWarp.cpp', 
'libsrc/Sponge.cpp', 
'libsrc/SSRigger.cpp', 
'libsrc/StateValue.cpp', 
'libsrc/SVDBond.cpp', 
'libsrc/TextManager.cpp', 
'libsrc/Timer.cpp', 
'libsrc/Twist.cpp', 
'libsrc/VBondReader.cpp', 
'libsrc/vec3.cpp', 
'libsrc/WaterNetwork.cpp', 
'libsrc/WeightedMap.cpp', 
'libsrc/Whack.cpp')

lib_vag = library('vag', gen_src, base_files,
dependencies: [png_dep, fftw_dep, dep_vgeom, dep_ica, dep_ccp4, boost_dep], cpp_args: ['-w', '-std=c++11'], install: true)
inc_vag = include_directories('libsrc')
dep_vag = declare_dependency(link_with: lib_vag, include_directories: inc_vag)

# Command line interface

executable('vagabond', 'libsrc/main.cpp', dependencies: [dep_vag, boost_dep], cpp_args: ['-std=c++11'], install: true)
#executable('vagabond', 'libsrc/main.cpp', base_files, dependencies: [boost_dep, dep_ica, png_dep, fftw_dep, dep_vgeom, dep_ccp4], cpp_args: ['-std=c++11', '-pg'], link_args: ['-pg'], install: true)

# Graphical user interface
# (including all the Qt and OpenGL stuff)

# Pre-processing by Qt5

if qt5_dep.found() and dep_gl.found()
	moc_files = qt5.preprocess(moc_headers : ['libgui/qtgui/CrystalExplorer.h',
'libgui/qtgui/Dialogue.h',
'libgui/qtgui/DropDown.h',
'libgui/qtgui/ErroneousZone.h',
'libgui/qtgui/MoleculeExplorer.h',
'libgui/qtgui/MonomerExplorer.h',
'libgui/qtgui/ResButton.h',
'libgui/qtgui/SequenceView.h',
'libgui/qtgui/SetterEdit.h',
'libgui/qtgui/StartScreen.h',
'libgui/qtgui/VagWindow.h'],
			moc_extra_arguments: ['-DMAKES_MY_MOC_HEADER_COMPILE'])

# Non-Qt5 GUI elements (mostly OpenGL stuff)
# Uses Qt5-supplied OpenGL

executable('vagabond-gui', base_files,
'libgui/Atoms2GL.cpp', 
'libgui/Bonds2GL.cpp', 
'libgui/Connect2GL.cpp', 
'libgui/Density2GL.cpp',
'libgui/GLKeeper.cpp',
'libgui/GLObject.cpp',
'libgui/Multi2GL.cpp', 
'libgui/shader.cpp',
'libgui/Selected2GL.cpp',
'libgui/Vagabond2GL.cpp', 
'libgui/WarpGL.cpp', 
'libgui/qtgui/gui_main.cpp',
'libgui/qtgui/CrystalExplorer.cpp',
'libgui/qtgui/Dialogue.cpp',
'libgui/qtgui/DropDown.cpp',
'libgui/qtgui/ErroneousZone.cpp',
'libgui/qtgui/gui_main.cpp',
'libgui/qtgui/InstructionThread.cpp',
'libgui/qtgui/MoleculeExplorer.cpp',
'libgui/qtgui/MonomerExplorer.cpp',
'libgui/qtgui/ResButton.cpp',
'libgui/qtgui/SequenceView.cpp',
'libgui/qtgui/SetterEdit.cpp',
'libgui/qtgui/StartScreen.cpp',
'libgui/qtgui/VagabondGLWidget.cpp',
'libgui/qtgui/VagWindow.cpp', 
moc_files, dependencies: [qt5_dep, dep_ica, dep_ccp4, dep_vgeom, png_dep, fftw_dep, boost_dep], install: true, cpp_args: ['-w',
'-std=c++11'])
else
	message('QT5 not found, so not building vagabond-gui.')
endif

# Pre-processing by Qt5

if qt5_dep.found() and dep_gl.found()
	moc_cluster4x = qt5.preprocess(moc_headers : [
'c4xsrc/Averager.h', 
'c4xsrc/AxisScroll.h', 
'c4xsrc/HKLView.h', 
'c4xsrc/GLPoint.h', 
'c4xsrc/Screen.h', 
'c4xsrc/SelectionWindow.h', 
'c4xsrc/KeeperGL.h', 
'c4xsrc/ClusterList.h', 
'c4xsrc/CorrelLabel.h', 
],
		moc_extra_arguments: ['-DMAKES_MY_MOC_HEADER_COMPILE'])


executable('force_down', gen_src, 
'force_down/main.cpp', 
'force_down/Wilson.cpp', 
cpp_args: ['-std=c++11'], dependencies : [ dep_vag, png_dep, dep_ica, dep_ccp4  ], install: true)

executable('cluster4x', gen_src, moc_cluster4x,
'c4xsrc/main.cpp', 
'c4xsrc/Averager.cpp', 
'c4xsrc/AxisScroll.cpp', 
'c4xsrc/ClusterList.cpp', 
'c4xsrc/CorrelLabel.cpp', 
'c4xsrc/MatrixView.cpp', 
'c4xsrc/MtzFile.cpp', 
'c4xsrc/MtzFFT.cpp', 
'c4xsrc/Screen.cpp', 
'c4xsrc/KeeperGL.cpp', 
'c4xsrc/HKLView.cpp', 
'c4xsrc/GLAxis.cpp', 
'c4xsrc/GLPoint.cpp', 
'c4xsrc/GLObject.cpp', 
'c4xsrc/SelectionWindow.cpp', 
cpp_args: ['-std=c++11'], dependencies : [ qt5_dep, dep_gl, dep_vag, png_dep, dep_ica ], install: true)
else
	message('QT5 not found, so not building cluster4x.')
endif


pkg = import('pkgconfig')
libs = [lib_vag, lib_ica]
h = ['libsrc'] # subdirectories of ${prefix}/${includedir} to add to header path
pkg.generate(libraries : libs,
             subdirs : h,
             version : '1.0',
             name : 'libvag',
             filebase : 'simple',
             description : 'Library for vagabond')

