set path=../c4xsrc/**,

command! Tags !ctags -R .
command! Ninja :wa|!ninja -C ../build/current
command! Dinja :wa|!ninja -C ../build/debug



