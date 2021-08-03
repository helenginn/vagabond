here = join_paths(meson.source_root(), 'install-scripts/cluster4x.app')

install_data('./res/Info.plist', 
		install_dir: join_paths(here, 'Contents'))
install_data('./res/HelloQt.icns', 
		install_dir: join_paths(here, 'Contents/Resources'))
install_data('./res/qt.conf', 
		install_dir: join_paths(here, 'Contents/Resources'))

meson.add_install_script('install_script.sh')
