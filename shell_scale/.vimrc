set path=../force_down/**,

command! Tags !ctags -R .
command! Ninja :wa|!ninja -C ../build/current



