project('vagabond', 'cpp', 'c')
qt5 = import('qt5')

qt5_dep = dependency('qt5', modules: ['Core', 'Gui', 'Widgets', 'Sql'], required : false)
dep_gl = dependency('gl', required : false)

png_dep = dependency('libpng')
boost_dep = dependency('boost')
fftw_dep = dependency('fftw3f')
maria_dep = dependency('libmariadb', required: false)

# Geometry and/or rotamer lookup madness

lib_vgeom = library('vgeom', 
'libinfo/GeomTable.cpp', 
'libinfo/RotamerTable.cpp', 
'libinfo/GeomVariant.cpp', 
cpp_args: ['-w'], install: true)
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
'libsrc/Angler.cpp', 
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
'libsrc/DistanceMatrix.cpp', 
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
'libsrc/SymMate.cpp', 
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

executable('force_down', gen_src, 
'force_down/main.cpp', 
'force_down/Wilson.cpp', 
cpp_args: ['-std=c++11'], dependencies : [ boost_dep, dep_vag, png_dep, dep_ica, dep_ccp4  ], install: true)


executable('shell_scale', gen_src, 
'shell_scale/main.cpp', 
cpp_args: ['-std=c++11'], dependencies : [ boost_dep, dep_vag, png_dep, dep_ica, dep_ccp4  ], install: true)

if qt5_dep.found() and dep_gl.found()

headers = [
'c4xsrc/Group.h', 
'c4xsrc/AxisScroll.h', 
'c4xsrc/CAlphaView.h', 
'c4xsrc/HKLView.h', 
'c4xsrc/FolderInput.h', 
'c4xsrc/Plot3D.h', 
'c4xsrc/PlotView.h', 
'c4xsrc/Screen.h', 
'c4xsrc/SelectionWindow.h', 
'c4xsrc/SelectMatrix.h', 
'c4xsrc/KeeperGL.h', 
'c4xsrc/ClusterList.h', 
'c4xsrc/CorrelLabel.h'
]

install_headers([
'libsrc/Absolute.h',
'libsrc/Anchor.h',
'libsrc/Angler.h',
'libsrc/Anisotropicator.h',
'libsrc/Any.h',
'libsrc/Atom.h',
'libsrc/AtomGroup.h',
'libsrc/Backbone.h',
'libsrc/Balance.h',
'libsrc/BaseParser.h',
'libsrc/Bond.h',
'libsrc/BondGroup.h',
'libsrc/Bucket.h',
'libsrc/BucketBulkSolvent.h',
'libsrc/BucketFeatures.h',
'libsrc/BucketPerStrand.h',
'libsrc/CSV.h',
'libsrc/Chelate.h',
'libsrc/Converter.h',
'libsrc/Crystal.h',
'libsrc/Dataset.h',
'libsrc/Diffraction.h',
'libsrc/DiffractionMTZ.h',
'libsrc/DistanceMatrix.h',
'libsrc/Element.h',
'libsrc/ExplicitModel.h',
'libsrc/FFT.h',
'libsrc/Fibonacci.h',
'libsrc/FileReader.h',
'libsrc/FlexGlobal.h',
'libsrc/FlexLocal.h',
'libsrc/GhostBond.h',
'libsrc/Hydrogenator.h',
'libsrc/KeyPoints.h',
'libsrc/Knotter.h',
'libsrc/MapScoreWorkspace.h',
'libsrc/Model.h',
'libsrc/Molecule.h',
'libsrc/Monomer.h',
'libsrc/Motion.h',
'libsrc/Notifiable.h',
'libsrc/Novalent.h',
'libsrc/Object.h',
'libsrc/Options.h',
'libsrc/PDBReader.h',
'libsrc/PNGFile.h',
'libsrc/Param.h',
'libsrc/ParamBand.h',
'libsrc/Parser.h',
'libsrc/ParserTypes.h',
'libsrc/PartialStructure.h',
'libsrc/Polymer.h',
'libsrc/Protocol.h',
'libsrc/Quat4Refine.h',
'libsrc/RefineMat3x3.h',
'libsrc/RefinementGridSearch.h',
'libsrc/RefinementLBFGS.h',
'libsrc/RefinementList.h',
'libsrc/RefinementNelderMead.h',
'libsrc/RefinementStepSearch.h',
'libsrc/RefinementStrategy.h',
'libsrc/Refitter.h',
'libsrc/SSRigger.h',
'libsrc/SVDBond.h',
'libsrc/Sampler.h',
'libsrc/Shouter.h',
'libsrc/Sidechain.h',
'libsrc/SpaceWarp.h',
'libsrc/Sponge.h',
'libsrc/StateValue.h',
'libsrc/SymMate.h',
'libsrc/TextManager.h',
'libsrc/Timer.h',
'libsrc/Twist.h',
'libsrc/VBondReader.h',
'libsrc/WaterNetwork.h',
'libsrc/WeightedMap.h',
'libsrc/Whack.h',
'libsrc/charmanip.h',
'libsrc/font.h',
'libsrc/mat3x3.h',
'libsrc/mat4x4.h',
'libsrc/maths.h',
'libsrc/polyfit.h',
'libsrc/quat4.h',
'libsrc/shared_ptrs.h',
'libsrc/vagcout.h',
'libsrc/vec3.h'],
 subdir : 'libsrc')

install_headers([
'c4xsrc/AveCAlpha.h', 
'c4xsrc/AveDiffraction.h', 
'c4xsrc/AveUnitCell.h', 
'c4xsrc/Average.h', 
'c4xsrc/AxisScroll.h', 
'c4xsrc/CAlphaView.h', 
'c4xsrc/ClusterList.h', 
'c4xsrc/CorrelLabel.h', 
'c4xsrc/DatasetPath.h',
'c4xsrc/AveCSV.h', 
'c4xsrc/FolderInput.h', 
'c4xsrc/GLAxis.h', 
'c4xsrc/GLPoint.h', 
'c4xsrc/GLObject.h', 
'c4xsrc/Group.h', 
'c4xsrc/HKLView.h', 
'c4xsrc/KeeperGL.h', 
'c4xsrc/MatrixView.h', 
'c4xsrc/MtzFile.h', 
'c4xsrc/MtzFFT.h', 
'c4xsrc/MtzFFTPtr.h', 
'c4xsrc/Plot3D.h', 
'c4xsrc/PlotR.h',
'c4xsrc/PlotView.h', 
'c4xsrc/QuickAtoms.h',
'c4xsrc/Output.h',
'c4xsrc/Screen.h', 
'c4xsrc/SelectionWindow.h',
'c4xsrc/SelectMatrix.h', 
'c4xsrc/UCPlot.h'],
 subdir : 'c4xsrc')

install_headers([
'liblbfgs/lbfgs.h'], 
 subdir : 'liblbfgs')

install_headers([
'c4xsrc/shaders/Blob_vsh.h', 
'c4xsrc/shaders/Blob_fsh.h'],
 subdir : 'c4xsrc/shaders')

install_headers([
'libccp4/cmtzlib.h', 
'libccp4/csymlib.h',
'libccp4/ccp4_spg.h',
'libccp4/ccp4_general.h'],
 subdir : 'libccp4')

c4xsrc = [
'c4xsrc/AveCAlpha.cpp', 
'c4xsrc/AveDiffraction.cpp', 
'c4xsrc/AveUnitCell.cpp', 
'c4xsrc/Average.cpp', 
'c4xsrc/AxisScroll.cpp', 
'c4xsrc/CAlphaView.cpp', 
'c4xsrc/ClusterList.cpp', 
'c4xsrc/CorrelLabel.cpp', 
'c4xsrc/AveCSV.cpp', 
'c4xsrc/FolderInput.cpp', 
'c4xsrc/GLAxis.cpp', 
'c4xsrc/GLPoint.cpp', 
'c4xsrc/GLObject.cpp', 
'c4xsrc/Group.cpp', 
'c4xsrc/HKLView.cpp', 
'c4xsrc/KeeperGL.cpp', 
'c4xsrc/MatrixView.cpp', 
'c4xsrc/MtzFile.cpp', 
'c4xsrc/MtzFFT.cpp', 
'c4xsrc/Plot3D.cpp', 
'c4xsrc/PlotR.cpp',
'c4xsrc/PlotView.cpp', 
'c4xsrc/QuickAtoms.cpp',
'c4xsrc/Output.cpp',
'c4xsrc/Screen.cpp', 
'c4xsrc/SelectionWindow.cpp',
'c4xsrc/SelectMatrix.cpp', 
'c4xsrc/UCPlot.cpp'
]


c4xdeps = [
qt5_dep,
dep_gl,
dep_vag,
dep_ica,
png_dep
]

if (maria_dep.found())
c4xsrc += ['c4xsrc/ClusterListSQL.cpp',
'c4xsrc/SQLCredentials.cpp', 
'c4xsrc/SQLInput.cpp', 
'c4xsrc/OutputSQL.cpp', 
'c4xsrc/Query.cpp' ]

headers += [
'c4xsrc/SQLCredentials.h', 
'c4xsrc/SQLInput.h'
]

c4xdeps += maria_dep

else

c4xsrc += 'c4xsrc/ClusterListNoSQL.cpp'
c4xsrc += 'c4xsrc/OutputNoSQL.cpp'
endif

moc_cluster4x = qt5.preprocess(moc_headers : headers,
		moc_extra_arguments: ['-DMAKES_MY_MOC_HEADER_COMPILE'])

lib_c4x = library('cluster4x', c4xsrc, moc_cluster4x, 
dependencies: c4xdeps, cpp_args: ['-Wall', '-std=c++11'], install: true)
inc_c4x = include_directories('c4xsrc')
dep_c4x = declare_dependency(link_with: lib_c4x, include_directories: inc_c4x)


executable('cluster4x', ['c4xsrc/main.cpp'], 
cpp_args: ['-std=c++11', '-Wall'], dependencies : [c4xdeps, dep_c4x], install: true)
else
	message('QT5 not found, so not building cluster4x.')
endif


pkg = import('pkgconfig')
h = ['libsrc', 'c4xsrc'] # subdirectories of ${prefix}/${includedir} to add to header path
pkg.generate(libraries : [lib_vag],
             subdirs : ['libsrc'],
             version : '1.0',
             name : 'libvag',
             filebase : 'libvag',
             description : 'Library for vagabond')
pkg.generate(libraries : [lib_ica],
             subdirs : ['libsrc'],
             version : '1.0',
             name : 'libica',
             filebase : 'libica',
             description : 'Library for ica and lbfgs')
pkg.generate(libraries : [lib_c4x],
             subdirs : ['c4xsrc'],
             version : '1.0',
             name : 'libcluster4x',
             filebase : 'libcluster4x',
             description : 'Library for cluster4x')


