set path+=libsrc/**,libgui/**,libinfo/**

command! Tags !ctags -R libgui/* libsrc/*
command! Ninja :!ninja -C build/current
command! Make1 !cd libgui/qtgui; make;

command! Doxy !doxygen Doxyfile

command Indent normal 0ggggVG=

" Unused at present.
":let mapleader = ","


