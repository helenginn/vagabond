######################################################################
# Automatically generated by qmake (3.1) Sun Dec 3 19:03:56 2017
######################################################################

TEMPLATE = app
TARGET = vagabond-gui.app
INCLUDEPATH += /usr/local/include /usr/local/opt/qt/include ../../
LIBS += -L../../libs/.libs/ -L../../libfftw/.libs -L/usr/local/lib -L/usr/lib64
LIBS += -lpng -lfftw3f -lica -lccp4c
QMAKE_LFLAGS += -Wl,-rpath,/usr/local/lib
QMAKE_CXXFLAGS += -O3 -g
QT += widgets
QT += opengl

# The following define makes your compiler warn you if you use any
# feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

# Input
HEADERS += ../../libsrc/Absolute.h ../../libsrc/Anchor.h ../../libsrc/Anisotropicator.h ../../libsrc/Atom.h ../../libsrc/AtomGroup.h ../../libsrc/Backbone.h ../../libsrc/Bond.h ../../libsrc/Bucket.h ../../libsrc/BucketUniform.h ../../libsrc/CSV.h ../../libsrc/Crystal.h ../../libsrc/Dataset.h ../../libsrc/Diffraction.h ../../libsrc/DiffractionMTZ.h ../../libsrc/Distributor.h ../../libsrc/Element.h ../../libsrc/FileReader.h ../../libsrc/FlexGlobal.h ../../libsrc/Kabsch.h ../../libsrc/Knotter.h ../../libsrc/LocalCC.h ../../libsrc/Model.h ../../libsrc/Molecule.h ../../libsrc/Monomer.h ../../libsrc/Object.h ../../libsrc/Options.h ../../libsrc/PDBReader.h ../../libsrc/PNGFile.h ../../libsrc/Polymer.h ../../libsrc/RefinementGridSearch.h ../../libsrc/RefinementNelderMead.h ../../libsrc/RefinementStepSearch.h ../../libsrc/RefinementStrategy.h ../../libsrc/Sampler.h ../../libsrc/Sandbox.h ../../libsrc/Shouter.h ../../libsrc/Sidechain.h ../../libsrc/TextManager.h ../../libsrc/fftw3d.h ../../libsrc/font.h ../../libsrc/mat3x3.h ../../libsrc/mat4x4.h ../../libsrc/maths.h ../../libsrc/shared_ptrs.h ../../libsrc/vec3.h VagWindow.h VagabondGLWidget.h ../../libsrc/Notifiable.h InstructionThread.h Dialogue.h MoleculeExplorer.h SequenceView.h ResButton.h MonomerExplorer.h SetterEdit.h ../../libsrc/Parser.h ../../libsrc/VBondReader.h ../../libinfo/GeomTable.h ../../libinfo/RotamerTable.h ../../libsrc/SSRigger.h ../Shaders/Shader_vsh.h ../Shaders/Shader_fsh.h
SOURCES += ../../libsrc/Absolute.cpp ../../libsrc/Anchor.cpp ../../libsrc/Anisotropicator.cpp ../../libsrc/Atom.cpp ../../libsrc/AtomGroup.cpp ../../libsrc/Backbone.cpp ../../libsrc/Bond.cpp ../../libsrc/Bucket.cpp ../../libsrc/BucketUniform.cpp ../../libsrc/CSV.cpp ../../libsrc/Crystal.cpp ../../libsrc/Diffraction.cpp ../../libsrc/DiffractionMTZ.cpp ../../libsrc/Distributor.cpp ../../libsrc/Element.cpp ../../libsrc/FileReader.cpp ../../libsrc/FlexGlobal.cpp ../../libsrc/Kabsch.cpp ../../libsrc/Knotter.cpp ../../libsrc/LocalCC.cpp ../../libsrc/Model.cpp ../../libsrc/Molecule.cpp ../../libsrc/Monomer.cpp ../../libsrc/Options.cpp ../../libsrc/PDBReader.cpp ../../libsrc/PNGFile.cpp ../../libsrc/Polymer.cpp ../../libsrc/RefinementGridSearch.cpp ../../libsrc/RefinementNelderMead.cpp ../../libsrc/RefinementStepSearch.cpp ../../libsrc/RefinementStrategy.cpp ../../libsrc/Sampler.cpp ../../libsrc/Sandbox.cpp ../../libsrc/Shouter.cpp ../../libsrc/Sidechain.cpp ../../libsrc/TextManager.cpp ../../libsrc/fftw3d.cpp ../../libsrc/gui_main.cpp ../../libsrc/mat3x3.cpp ../../libsrc/mat4x4.cpp ../../libsrc/maths.cpp ../../libsrc/vec3.cpp VagWindow.cpp VagabondGLWidget.cpp ../GLKeeper.cpp ../GLObject.cpp ../Vagabond2GL.cpp ../shader.cpp InstructionThread.cpp Dialogue.cpp MoleculeExplorer.cpp SequenceView.cpp ResButton.cpp MonomerExplorer.cpp SetterEdit.cpp ../../libsrc/Parser.cpp ../../libsrc/VBondReader.cpp ../../libinfo/GeomTable.cpp ../../libinfo/RotamerTable.cpp ../../libsrc/SSRigger.cpp
