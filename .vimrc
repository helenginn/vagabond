set path+=libsrc/**,libgui/**,libinfo/**

command! MakeTags !ctags -R libgui/* libsrc/*
command! Make !cd libgui/qtgui; make -j;
