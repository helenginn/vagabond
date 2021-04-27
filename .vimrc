set path=libsrc/**,libgui/**,libinfo/**,doc/,c4xsrc/**,subprojects/helen3d/h3dsrc/**,subprojects/helencore/hcsrc/**,force_down/**

command! Tags !ctags -R libgui/* libsrc/*
command! Ninja :wa|!ninja -C build/current
command! Make1 !cd libgui/qtgui; make;

command! Doxy !doxygen Doxyfile



