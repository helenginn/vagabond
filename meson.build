project('vagabond', 'cpp', 'c')
qt5 = import('qt5')
cxx = meson.get_compiler('cpp')
arg_list = []

if (host_machine.system() == 'darwin')
  arg_list += '-mmacosx-version-min=10.14'
endif

add_global_arguments(arg_list, language : ['cpp', 'c'])
add_global_arguments('-std=c++11', language : ['cpp'])


qt5_dep = dependency('qt5', modules: ['Core', 'Gui', 'Widgets', 'Sql'], required : false)
dep_gl = dependency('gl', required : false)

png_dep = dependency('libpng')
boost_dep = dependency('boost')
fftw_dep = dependency('fftw3f')
thread_fft_dep = dependency('fftw3f')
thread_dep = dependency('threads')
maria_dep = dependency('libmariadb', required: false)
crystfel_dep = dependency('crystfel', required: false)

test_fftwf_threads = '''
#include <fftw3.h>
int main() {
    fftwf_init_threads();
    return 0;
}
'''
if not cxx.links(test_fftwf_threads, dependencies: fftw_dep)
  thread_fft_dep = cxx.find_library('fftw3f_threads')
  message('find_library')
endif

sp_helen3d = subproject('helen3d')
helen3d_dep = dependency('', required: false)

if (sp_helen3d.get_variable('helen3d_scs'))
	helen3d_dep = sp_helen3d.get_variable('helen3d_subdep')
endif

helencore_dep = sp_helen3d.get_variable('helencore_dep')

# Geometry and/or rotamer lookup madness

lib_vgeom = library('vgeom', 
'libinfo/GeomTable.cpp', 
'libinfo/RotamerTable.cpp', 
'libinfo/GeomVariant.cpp', 
install_rpath: '$ORIGIN/.',
cpp_args: arg_list, link_args: arg_list, install: true)
inc_vgeom = include_directories('libinfo')
dep_vgeom = declare_dependency(link_with: lib_vgeom, include_directories: inc_vgeom)

# numerical recipes, single value decomposition

lib_ica = library('ica', 'libica/matrix.cpp', 'libica/svdcmp.cpp', 'liblbfgs/lbfgs.c', 
install_rpath: '$ORIGIN/.',
cpp_args: arg_list, link_args: arg_list, install: true)
inc_ica = include_directories('libica')
dep_ica = declare_dependency(link_with: lib_ica, include_directories: inc_ica)

# cmtzlib / csymlib: SYMINFO-variable-independent distributions needs separate compiling
# May or may not need the maths library on various systems

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)

lib_ccp4 = library('ccp4', 'libccp4/ccp4_array.c', 'libccp4/ccp4_general.c', 'libccp4/ccp4_parser.c', 'libccp4/ccp4_program.c', 'libccp4/ccp4_unitcell.c', 'libccp4/cmap_accessor.c', 'libccp4/cmap_close.c', 'libccp4/cmap_data.c', 'libccp4/cmap_header.c', 'libccp4/cmap_labels.c', 'libccp4/cmap_open.c', 'libccp4/cmap_skew.c', 'libccp4/cmap_stats.c', 'libccp4/cmap_symop.c', 'libccp4/cmtzlib.c', 'libccp4/csymlib.c', 'libccp4/cvecmat.c', 'libccp4/library_err.c', 'libccp4/library_file.c', 'libccp4/library_utils.c', 'libccp4/pack_c.c', dependencies: m_dep, 
cpp_args: arg_list, link_args: arg_list, install: true)
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
'libsrc/Blob.cpp', 
'libsrc/Bond.cpp', 
'libsrc/BondGroup.cpp', 
'libsrc/BucketBulkSolvent.cpp', 
'libsrc/BucketPerStrand.cpp', 
'libsrc/Bucket.cpp', 
'libsrc/Chelate.cpp', 
'libsrc/CompareFFT.cpp', 
'libsrc/Crystal.cpp', 
'libsrc/CSV.cpp', 
'libsrc/Diffraction.cpp', 
'libsrc/DiffractionMTZ.cpp', 
'libsrc/DistanceMatrix.cpp', 
'libsrc/Element.cpp', 
'libsrc/ExplicitModel.cpp', 
'libsrc/FFT.cpp', 
'libsrc/FlexLocal.cpp', 
'libsrc/GhostBond.cpp',
'libsrc/Hydrogenator.cpp', 
'libsrc/KeyPoints.cpp', 
'libsrc/Knotter.cpp', 
'libsrc/Model.cpp', 
'libsrc/Molecule.cpp', 
'libsrc/Motion.cpp', 
'libsrc/Monomer.cpp', 
'libsrc/Multistate.cpp', 
'libsrc/Novalent.cpp', 
'libsrc/Options.cpp', 
'libsrc/ParamBand.cpp', 
'libsrc/Parser.cpp', 
'libsrc/PartialStructure.cpp', 
'libsrc/PDBReader.cpp', 
'libsrc/PNGFile.cpp', 
'libsrc/Polymer.cpp', 
'libsrc/Refitter.cpp', 
'libsrc/Sampler.cpp', 
'libsrc/Scaler.cpp', 
'libsrc/Shouter.cpp', 
'libsrc/Sidechain.cpp', 
'libsrc/Sponge.cpp', 
'libsrc/SSRigger.cpp', 
'libsrc/StateValue.cpp', 
'libsrc/SymAtom.cpp', 
'libsrc/SymMate.cpp', 
'libsrc/SVDBond.cpp', 
'libsrc/TextManager.cpp', 
'libsrc/Twist.cpp', 
'libsrc/VBondReader.cpp', 
'libsrc/vagmaths.cpp', 
'libsrc/WaterNetwork.cpp', 
'libsrc/WeightedMap.cpp', 
'libsrc/Whack.cpp')

lib_vag = library('vag', gen_src, base_files,
include_directories: 'subprojects/helencore',
install_rpath: '$ORIGIN/.',
dependencies: [png_dep, fftw_dep, thread_dep, thread_fft_dep, dep_vgeom, dep_ica, dep_ccp4, boost_dep,
helencore_dep], cpp_args: ['-w', '-std=c++11'], link_args: arg_list, install: true)
inc_vag = include_directories('libsrc')
dep_vag = declare_dependency(link_with: lib_vag, include_directories: inc_vag)

# Command line interface

executable('vagabond', 'libsrc/main.cpp', include_directories: 'subprojects/helencore', dependencies: [dep_vag, boost_dep], cpp_args: ['-std=c++11'], install: true)
#executable('vagabond', 'libsrc/main.cpp', base_files, dependencies: [boost_dep, dep_ica, png_dep, fftw_dep, dep_vgeom, dep_ccp4], cpp_args: ['-std=c++11', '-pg'], link_args: ['-pg'], install: true)

# Graphical user interface
# (including all the Qt and OpenGL stuff)

# Pre-processing by Qt5

if qt5_dep.found() and dep_gl.found()
	moc_files = qt5.preprocess(moc_headers : ['libgui/qtgui/CrystalExplorer.h',
'libgui/qtgui/Dialogue.h',
'libgui/qtgui/DropDown.h',
'libgui/qtgui/MoleculeExplorer.h',
'libgui/qtgui/MonomerExplorer.h',
'libgui/qtgui/ResButton.h',
'libgui/qtgui/SequenceView.h',
'libgui/qtgui/SetterEdit.h',
'libgui/qtgui/StartScreen.h',
'libgui/qtgui/VagabondGLWidget.h',
'libgui/qtgui/VagWindow.h'],
			moc_extra_arguments: ['-DMAKES_MY_MOC_HEADER_COMPILE'])

# Non-Qt5 GUI elements (mostly OpenGL stuff)
# Uses Qt5-supplied OpenGL

subincludes = ['subprojects/helen3d', 'subprojects/helencore']

files = [
'libgui/Atoms2GL.cpp', 
'libgui/BlobMesh.cpp', 
'libgui/Bonds2GL.cpp', 
'libgui/Connect2GL.cpp', 
'libgui/Density2GL.cpp',
'libgui/Multi2GL.cpp', 
'libgui/Selected2GL.cpp',
'libgui/Vagabond2GL.cpp', 
'libgui/qtgui/gui_main.cpp',
'libgui/qtgui/CrystalExplorer.cpp',
'libgui/qtgui/Dialogue.cpp',
'libgui/qtgui/DropDown.cpp',
'libgui/qtgui/gui_main.cpp',
'libgui/qtgui/InstructionThread.cpp',
'libgui/qtgui/MoleculeExplorer.cpp',
'libgui/qtgui/MonomerExplorer.cpp',
'libgui/qtgui/ResButton.cpp',
'libgui/qtgui/SequenceView.cpp',
'libgui/qtgui/SetterEdit.cpp',
'libgui/qtgui/StartScreen.cpp',
'libgui/qtgui/VagabondGLWidget.cpp',
'libgui/qtgui/VagWindow.cpp'
]

executable('vagabond-gui', base_files, files,
moc_files, dependencies: [qt5_dep, dep_ica, helen3d_dep, helencore_dep,
dep_ccp4, dep_vgeom, png_dep, fftw_dep, thread_dep, thread_fft_dep, boost_dep], 
include_directories: subincludes,
install: true, cpp_args: ['-w', '-std=c++11'])

lib_vgui = library('vgui', base_files, files, moc_files,
dependencies: [qt5_dep, dep_ica, helen3d_dep, helencore_dep,
dep_ccp4, dep_vgeom, png_dep, fftw_dep, thread_dep, thread_fft_dep, boost_dep], 
install_rpath: '$ORIGIN/.',
include_directories: subincludes,
cpp_args: ['-w', '-std=c++11'], install: true)
else
	message('QT5 not found, so not building vagabond-gui.')
endif


executable('force_down', gen_src, 
'force_down/main.cpp', 
'force_down/Wilson.cpp', 
cpp_args: ['-std=c++11'], 
include_directories: 'subprojects/helencore',
dependencies : [ helencore_dep, boost_dep, dep_vag, png_dep, dep_ica, dep_ccp4  ], install: true)


executable('shell_scale', gen_src, 
'shell_scale/main.cpp', 
cpp_args: ['-std=c++11'],
include_directories: 'subprojects/helencore',
dependencies : [ helencore_dep, boost_dep, dep_vag, png_dep, dep_ica, dep_ccp4  ], install: true)

if qt5_dep.found() and dep_gl.found()

headers = [
'c4xsrc/Arbitrary.h', 
'c4xsrc/Group.h', 
'c4xsrc/AxisScroll.h', 
'c4xsrc/CAlphaView.h', 
'c4xsrc/ColumnView.h', 
'c4xsrc/HKLView.h', 
'c4xsrc/FolderInput.h', 
'c4xsrc/Input.h', 
'c4xsrc/ClusterPlot.h', 
'c4xsrc/DisplaySettings.h', 
'c4xsrc/PlotView.h', 
'c4xsrc/Screen.h', 
'c4xsrc/SelectionWindow.h', 
'c4xsrc/SelectMatrix.h', 
'c4xsrc/KeeperGL.h', 
'c4xsrc/ClusterList.h', 
'c4xsrc/CorrelLabel.h',
'c4xsrc/MtzOptions.h', 
'c4xsrc/AutoCluster.h', 
]

install_headers([
'libsrc/Absolute.h',
'libsrc/Anchor.h',
'libsrc/Angler.h',
'libsrc/Anisotropicator.h',
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
'libsrc/Crystal.h',
'libsrc/Dataset.h',
'libsrc/Diffraction.h',
'libsrc/DiffractionMTZ.h',
'libsrc/DistanceMatrix.h',
'libsrc/Element.h',
'libsrc/ExplicitModel.h',
'libsrc/FFT.h',
'libsrc/FlexLocal.h',
'libsrc/GhostBond.h',
'libsrc/Hydrogenator.h',
'libsrc/KeyPoints.h',
'libsrc/Knotter.h',
'libsrc/lbfgs.h',
'libsrc/MapScoreWorkspace.h',
'libsrc/Model.h',
'libsrc/Molecule.h',
'libsrc/Monomer.h',
'libsrc/Multistate.h',
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
'libsrc/Refitter.h',
'libsrc/SSRigger.h',
'libsrc/SVDBond.h',
'libsrc/Sampler.h',
'libsrc/Shouter.h',
'libsrc/Sidechain.h',
'libsrc/Sponge.h',
'libsrc/StateValue.h',
'libsrc/SymAtom.h',
'libsrc/SymMate.h',
'libsrc/TextManager.h',
'libsrc/Twist.h',
'libsrc/VBondReader.h',
'libsrc/WaterNetwork.h',
'libsrc/WeightedMap.h',
'libsrc/Whack.h',
'libsrc/charmanip.h',
'libsrc/font.h',
'libsrc/vagmaths.h',
'libsrc/polyfit.h',
'libsrc/quat4.h',
'libsrc/shared_ptrs.h',
'libsrc/vagcout.h'],
 subdir : 'libsrc')

install_headers([
'libgui/Density2GL.h',
],
 subdir : 'vaggui')

install_headers([
'libica/svdcmp.h',
'libica/matrix.h',
],
 subdir : 'libica')

install_headers([
'c4xsrc/AveCAlpha.h', 
'c4xsrc/AveDiffraction.h', 
'c4xsrc/AveUnitCell.h', 
'c4xsrc/Average.h', 
'c4xsrc/AxisScroll.h', 
'c4xsrc/CAlphaView.h', 
'c4xsrc/ClusterList.h', 
'c4xsrc/MtzOptions.h', 
'c4xsrc/CorrelLabel.h', 
'c4xsrc/DatasetPath.h',
'c4xsrc/AveCSV.h', 
'c4xsrc/AveVectors.h', 
'c4xsrc/FolderInput.h', 
'c4xsrc/ExportType.h', 
'c4xsrc/GLAxis.h', 
'c4xsrc/GLPoint.h', 
'c4xsrc/Group.h', 
'c4xsrc/HKLView.h', 
'c4xsrc/KeeperGL.h', 
'c4xsrc/MatrixView.h', 
'c4xsrc/MtzFile.h', 
'c4xsrc/MtzFFT.h', 
'c4xsrc/MtzFFTPtr.h', 
'c4xsrc/ClusterPlot.h', 
'c4xsrc/PlotR.h',
'c4xsrc/PlotView.h', 
'c4xsrc/QuickAtoms.h',
'c4xsrc/Output.h',
'c4xsrc/Screen.h', 
'c4xsrc/SelectionWindow.h',
'c4xsrc/SelectMatrix.h', 
'c4xsrc/Vertex.h', 
'c4xsrc/UCPlot.h'],
 subdir : 'c4xsrc')

install_headers([
'c4xsrc/shaders/Blob_vsh.h', 
'c4xsrc/shaders/Blob_fsh.h'],
 subdir : 'c4xsrc/shaders')

install_headers([
'liblbfgs/lbfgs.h'], 
 subdir : 'liblbfgs')

install_headers([
'c4xsrc/shaders/Blob_vsh.h', 
'c4xsrc/shaders/Blob_fsh.h'],
 subdir : 'c4xsrc/shaders')

install_headers([
'libccp4/ccp4_array.h',
'libccp4/ccp4_errno.h',
'libccp4/ccp4_file_err.h',
'libccp4/ccp4_fortran.h',
'libccp4/ccp4_general.h',
'libccp4/ccp4_parser.h',
'libccp4/ccp4_program.h',
'libccp4/ccp4_spg.h',
'libccp4/ccp4_sysdep.h',
'libccp4/ccp4_types.h',
'libccp4/ccp4_unitcell.h',
'libccp4/ccp4_utils.h',
'libccp4/ccp4_vars.h',
'libccp4/cmap_data.h',
'libccp4/cmap_errno.h',
'libccp4/cmap_header.h',
'libccp4/cmap_labels.h',
'libccp4/cmap_skew.h',
'libccp4/cmap_stats.h',
'libccp4/cmaplib.h',
'libccp4/cmaplib_f.h',
'libccp4/cmtzlib.h',
'libccp4/csymlib.h',
'libccp4/cvecmat.h',
'libccp4/library_file.h',
'libccp4/mtzdata.h',
'libccp4/overview.h',
'libccp4/pack_c.h',
'libccp4/syminfo.h',
'libccp4/w32mvs.h'],
 subdir : 'libccp4')

install_headers([
'libinfo/CentroidToPhase.h',
'libinfo/GeomTable.h',
'libinfo/GeomVariant.h',
'libinfo/GeomVariantDefs.h',
'libinfo/RotamerTable.h',
'libinfo/ScatterFactors.h'],
 subdir : 'libinfo')

c4xsrc = [
'c4xsrc/Arbitrary.cpp', 
'c4xsrc/AutoCluster.cpp', 
'c4xsrc/AveCAlpha.cpp', 
'c4xsrc/AveDiffraction.cpp', 
'c4xsrc/AveUnitCell.cpp', 
'c4xsrc/Average.cpp', 
'c4xsrc/AxisScroll.cpp', 
'c4xsrc/CAlphaView.cpp', 
'c4xsrc/ClusterList.cpp', 
'c4xsrc/ClusterPlot.cpp', 
'c4xsrc/ColumnView.cpp', 
'c4xsrc/ColumnChooser.cpp', 
'c4xsrc/CorrelLabel.cpp', 
'c4xsrc/DisplaySettings.cpp', 
'c4xsrc/AveCSV.cpp', 
'c4xsrc/AveVectors.cpp', 
'c4xsrc/FolderInput.cpp', 
'c4xsrc/GLAxis.cpp', 
'c4xsrc/GLPoint.cpp', 
'c4xsrc/Group.cpp', 
'c4xsrc/HKLView.cpp', 
'c4xsrc/Input.cpp', 
'c4xsrc/KeeperGL.cpp', 
'c4xsrc/MatrixView.cpp', 
'c4xsrc/MtzFile.cpp', 
'c4xsrc/MtzFFT.cpp', 
'c4xsrc/MtzOptions.cpp', 
'c4xsrc/MyDictator.cpp', 
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
dep_vag,
dep_ica,
png_dep,
dep_ccp4,
helencore_dep,
helen3d_dep
]

if (crystfel_dep.found())
c4xsrc += ['c4xsrc/CrystFELInput.cpp', ]
c4xsrc += ['c4xsrc/ClusterListCrystFEL.cpp', ]
c4xdeps += crystfel_dep

else
c4xsrc += ['c4xsrc/ClusterListNoCrystFEL.cpp', ]

endif

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
include_directories: subincludes,
install_rpath: '$ORIGIN/.',
dependencies: c4xdeps, cpp_args: ['-Wall', '-std=c++11'], link_args: arg_list, install: true)
inc_c4x = include_directories('c4xsrc')
dep_c4x = declare_dependency(link_with: lib_c4x, include_directories: [inc_c4x])

executable('cluster4x', ['c4xsrc/main.cpp'], 
cpp_args: ['-std=c++11', '-Wall'], 
include_directories : subincludes,
dependencies : [helencore_dep, c4xdeps, dep_c4x], link_args: arg_list, install: true)
else
	message('QT5 not found, so not building cluster4x.')
endif


pkg = import('pkgconfig')
h = ['libsrc', 'c4xsrc'] # subdirectories of ${prefix}/${includedir} to add to header path
pkg.generate([lib_vag],
             subdirs : ['libsrc'],
             version : '1.0',
             name : 'libvag',
             filebase : 'libvag',
             description : 'Library for vagabond')
if qt5_dep.found() and dep_gl.found()
pkg.generate([lib_vgui],
             subdirs : ['libgui'],
             version : '1.0',
             name : 'libvgui',
             filebase : 'libvgui',
             description : 'Library for vagabond GUI')
endif
pkg.generate([lib_ica],
             subdirs : ['libsrc'],
             version : '1.0',
             name : 'libica',
             filebase : 'libica',
             description : 'Library for ica and lbfgs')
pkg.generate([lib_vgeom],
             subdirs : ['libsrc'],
             version : '1.0',
             name : 'libvgeom',
             filebase : 'libvgeom',
             description : 'Library for vagabond residue info')
pkg.generate([lib_ccp4],
             subdirs : ['libsrc'],
             version : '1.0',
             name : 'libccp4',
             filebase : 'libccp4',
             description : 'Library for ccp4 (compiled using Vagabond)')

if (sp_helen3d.get_variable('helen3d_scs'))
pkg.generate([lib_c4x],
             subdirs : ['c4xsrc'],
             version : '1.0',
             name : 'libcluster4x',
             filebase : 'libcluster4x',
             description : 'Library for cluster4x')

endif

#subdir('install-scripts')
