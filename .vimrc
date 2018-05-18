set path+=libsrc/**,libgui/**,libinfo/**

command! MakeTags !ctags -R libgui/* libsrc/*
command! Make :make -j -C libgui/qtgui
command! Make1 !cd libgui/qtgui; make;

command! Doxy !doxygen Doxyfile

command Indent normal 0ggggVG=

:nnoremap <C-P> o<Esc>
