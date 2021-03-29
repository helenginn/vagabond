set path=../c4xsrc/**,../subprojects/helen3d/h3dsrc/**, ../subprojects/helencore/hcsrc/**

command! Tags !ctags -R .
command! Ninja :wa|!ninja -C ../build/current
command! Dinja :wa|!ninja -C ../build/debug



